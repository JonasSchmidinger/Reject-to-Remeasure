# In this script, we create Fig. S5
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(ggpubr)
library(prospectr)
library(ggpointdensity)
library(zoo)
library(ggh4x)
library(future.apply)
library(parallel)
library(metR)

plan(multisession, workers = 6) # Will be quite slow !!!


quantiles <- c(
  0.001,
  round(seq(0.005, 0.995, by = 0.005), 3),
  0.999
)


SOM_true<- read.csv(file = "Results/TabPFN/Outer/SOM_true.csv")
Clay_true<- read.csv(file = "Results/TabICL/Outer/Clay_true.csv")
K_Meh3_true<-  read.csv(file = "Results/TabPFN/Outer/K_Meh3_true.csv")
P_Meh3_true<-  read.csv(file = "Results/TabPFN/Outer/P_Meh3_true.csv")

SOM_pred<- read.csv(file = "Results/TabPFN/Outer/SOM_pred.csv")
Clay_pred<- read.csv(file = "Results/TabICL/Outer/Clay_pred.csv")
K_Meh3_pred<-  read.csv(file = "Results/TabPFN/Outer/K_Meh3_pred.csv")
P_Meh3_pred<-  read.csv(file = "Results/TabPFN/Outer/P_Meh3_pred.csv")

SOM_pred_quantiles<- read.csv(file = "Results/TabPFN/Outer/SOM_pred_quantiles.csv")
Clay_pred_quantiles<- read.csv(file = "Results/TabICL/Outer/Clay_pred_quantiles.csv")
K_Meh3_pred_quantiles<-  read.csv(file = "Results/TabPFN/Outer/K_Meh3_pred_quantiles.csv")
P_Meh3_pred_quantiles<-  read.csv(file = "Results/TabPFN/Outer/P_Meh3_pred_quantiles.csv")




compute_p_exceed_hybrid <- function(quantile_matrix,
                                    pred_values,
                                    a,
                                    b,
                                    quantile_probs) {

  n <- nrow(quantile_matrix)
  p_vals <- numeric(n)

  for (i in seq_len(n)) {

    q_vals <- as.numeric(quantile_matrix[i, ])

    ord <- order(q_vals)
    q_vals  <- q_vals[ord]
    q_probs <- quantile_probs[ord]

    cdf_fun <- stepfun(
      x = q_vals,
      y = c(0, q_probs),
      right = TRUE
    )

    y_hat <- pred_values[i]

    L <- (y_hat - a) / (1 + b)
    U <- (y_hat + a) / (1 - b)

    p_vals[i] <- cdf_fun(L) + (1 - cdf_fun(U))
  }

  data.frame(p_value = p_vals)
}


#----------------------------------------------------------------------------------------------------------------
# Clay

Clay_pred_values      <- as.numeric(Clay_pred$Clay_pred) * 10
Clay_true_values      <- as.numeric(Clay_true$Clay_target) * 10
Clay_quantile_matrix  <- as.matrix(Clay_pred_quantiles) * 10

Clay_error <- abs(Clay_pred_values - Clay_true_values)


delta0_values <-seq(-2, 152, length.out = 100)
b_values      <- seq(-0.02, 0.42, length.out = 100)

Clay_economic_df <- expand.grid(
  delta0 = delta0_values,
  b      = b_values
) %>%
  arrange(b, delta0)

Clay_economic_df$incorrect_acceptance_rate <- NA_real_

alpha_fixed <- 0.05

Clay_economic_df$incorrect_acceptance_rate <- future_sapply(
  seq_len(nrow(Clay_economic_df)),
  function(i) {

    a_i <- Clay_economic_df$delta0[i]
    b_i <- Clay_economic_df$b[i]

    p_vals <- compute_p_exceed_hybrid(
      quantile_matrix = Clay_quantile_matrix,
      pred_values     = Clay_pred_values,
      a               = a_i,
      b               = b_i,
      quantile_probs  = quantiles
    )

    delta_true <- a_i + b_i * Clay_true_values

    accepted <- p_vals$p_value < alpha_fixed
    wrong    <- Clay_error > delta_true

    if (sum(accepted) == 0) {
      0
    } else {
      mean(wrong[accepted])
    }

  },
  future.seed = TRUE
)

Clay_economic_df <- Clay_economic_df %>%
  mutate(
    incorrect_acceptance_percent = incorrect_acceptance_rate * 100
  )



Clay_economic_df <- Clay_economic_df %>%
  mutate(
    incorrect_acceptance_class = case_when(
      incorrect_acceptance_percent <= 2.5 ~ "0–2.5",
      incorrect_acceptance_percent <= 5   ~ "2.5–5",
      incorrect_acceptance_percent <= 7.5 ~ "5–7.5",
      incorrect_acceptance_percent <= 10  ~ "7.5–10",
      TRUE                                ~ ">10"
    )
  )

Clay_economic_df <- Clay_economic_df %>%
  mutate(
    incorrect_acceptance_class = factor(
      incorrect_acceptance_class,
      levels = c("0–2.5", "2.5–5", "5–7.5", "7.5–10", ">10")
    )
  )


Clay_economic_df$title <- "bold(Clay~(g~kg^{-1})*' - TabICL')"


Clay_Incorrect_Acceptance_Plot <- ggplot(
  Clay_economic_df,
  aes(
    x = delta0,
    y = b,
    fill = incorrect_acceptance_class
  )
) +

  geom_raster() +

  scale_fill_manual(
    values = c(
      "0–2.5" = "#1F4E79",
      "2.5–5" = "#B2E2E2",
      "5–7.5" = "#FEC44F",
      "7.5–10" = "#FB6A4A",
      ">10" = "#67000D"
    ),
    name = expression("Incorrect Acceptances Among Accepted Predictions (%)" * ", " * alpha == 0.05),
    drop = FALSE
  ) +

  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(
    xlim = c(0,150),
    ylim = c(0,0.4),
    expand = FALSE
  )+

  labs(
    x = expression(delta[0]),
    y = expression(italic(b))
  ) +

  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.title = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    legend.text = element_text(size = 18, colour = "black"),
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
    axis.line = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 18, colour = "black", face = "bold"),
    strip.background = element_rect(fill = "grey95", linewidth = 0.5)
  ) +
  facet_grid(. ~ title, labeller = label_parsed)

Clay_Incorrect_Acceptance_Plot


#----------------------------------------------------------------------------------------------------------------
# SOM


SOM_pred_values     <- as.numeric(SOM_pred$SOM_pred) * 10
SOM_true_values     <- as.numeric(SOM_true$SOM_target) * 10
SOM_quantile_matrix <- as.matrix(SOM_pred_quantiles) * 10

SOM_error <- abs(SOM_pred_values - SOM_true_values)

delta0_values <- seq(-0.05, 25.5, length.out = 100)
b_values      <- seq(-0.01, 0.42, length.out = 100)




SOM_economic_df <- expand.grid(
  delta0 = delta0_values,
  b      = b_values
) %>%
  arrange(b, delta0)

alpha_fixed <- 0.05

SOM_economic_df$incorrect_acceptance_rate <- future_sapply(
  seq_len(nrow(SOM_economic_df)),
  function(i) {

    a_i <- SOM_economic_df$delta0[i]
    b_i <- SOM_economic_df$b[i]

    p_vals <- compute_p_exceed_hybrid(
      quantile_matrix = SOM_quantile_matrix,
      pred_values     = SOM_pred_values,
      a               = a_i,
      b               = b_i,
      quantile_probs  = quantiles
    )

    delta_true <- a_i + b_i * SOM_true_values

    accepted <- p_vals$p_value < alpha_fixed
    wrong    <- SOM_error > delta_true

    if (sum(accepted) == 0) {
      0
    } else {
      mean(wrong[accepted])
    }
  },
  future.seed = TRUE
)

SOM_economic_df <- SOM_economic_df %>%
  mutate(
    incorrect_acceptance_percent = incorrect_acceptance_rate * 100
  )

SOM_economic_df

SOM_economic_df <- SOM_economic_df %>%
  mutate(
    incorrect_acceptance_class = case_when(
      incorrect_acceptance_percent <= 2.5 ~ "0–2.5",
      incorrect_acceptance_percent <= 5   ~ "2.5–5",
      incorrect_acceptance_percent <= 7.5 ~ "5–7.5",
      incorrect_acceptance_percent <= 10  ~ "7.5–10",
      TRUE                                ~ ">10"
    )
  )

SOM_economic_df <- SOM_economic_df %>%
  mutate(
    incorrect_acceptance_class = factor(
      incorrect_acceptance_class,
      levels = c("0–2.5", "2.5–5", "5–7.5", "7.5–10", ">10")
    )
  )

SOM_economic_df$title <- "bold(SOM~(g~kg^{-1})*' - TabPFN')"

SOM_Incorrect_Acceptance_Plot <- ggplot(
  SOM_economic_df,
  aes(
    x = delta0,
    y = b,
    fill = incorrect_acceptance_class
  )
) +
  geom_raster() +

  scale_fill_manual(
    values = c(
      "0–2.5" = "#1F4E79",
      "2.5–5" = "#B2E2E2",
      "5–7.5" = "#FEC44F",
      "7.5–10" = "#FB6A4A",
      ">10" = "#67000D"
    ),
    name = expression("Incorrect Acceptances Among Accepted Predictions (%)" * ", " * alpha == 0.05),
    drop = FALSE
  ) +

  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(
    xlim = c(0, 25),
    ylim = c(0, 0.4),
    expand = FALSE
  ) +

  labs(
    x = expression(delta[0]),
    y = expression(italic(b))
  ) +

  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    legend.position = c(-0.595, 1.18),
    legend.justification = c(0, 0),
    plot.margin = margin(t = 85, r = 10, b = 10, l = 10),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(size = 20, colour = "black"),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.text = element_text(size = 18, colour = "black", margin = margin(t = 3, r = 25)),
    panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
    axis.line = element_blank(),
    strip.text = element_text(size = 18, colour = "black", face = "bold"),
    strip.background = element_rect(fill = "grey95", linewidth = 0.5)
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      keywidth = unit(1.8, "cm"),
      keyheight = unit(0.8, "cm"),
      order = 1
    )
  ) +
  facet_grid(. ~ title, labeller = label_parsed)

SOM_Incorrect_Acceptance_Plot




#----------------------------------------------------------------------------------------------------------------
# K


K_Meh3_pred_values     <- as.numeric(K_Meh3_pred$K_Meh3_pred)
K_Meh3_true_values     <- as.numeric(K_Meh3_true$K_Meh3_target)
K_Meh3_quantile_matrix <- as.matrix(K_Meh3_pred_quantiles)

K_Meh3_error <- abs(K_Meh3_pred_values - K_Meh3_true_values)

delta0_values <- seq(-1, 81, length.out = 100)
b_values      <- seq(-0.01, 0.41, length.out = 100)

K_Meh3_economic_df <- expand.grid(
  delta0 = delta0_values,
  b      = b_values
) %>%
  arrange(b, delta0)

alpha_fixed <- 0.05

K_Meh3_economic_df$incorrect_acceptance_rate <- future_sapply(
  seq_len(nrow(K_Meh3_economic_df)),
  function(i) {

    a_i <- K_Meh3_economic_df$delta0[i]
    b_i <- K_Meh3_economic_df$b[i]

    p_vals <- compute_p_exceed_hybrid(
      quantile_matrix = K_Meh3_quantile_matrix,
      pred_values     = K_Meh3_pred_values,
      a               = a_i,
      b               = b_i,
      quantile_probs  = quantiles
    )

    delta_true <- a_i + b_i * K_Meh3_true_values

    accepted <- p_vals$p_value < alpha_fixed
    wrong    <- K_Meh3_error > delta_true

    if (sum(accepted) == 0) {
      0
    } else {
      mean(wrong[accepted])
    }

  },
  future.seed = TRUE
)

K_Meh3_economic_df <- K_Meh3_economic_df %>%
  mutate(
    incorrect_acceptance_percent = incorrect_acceptance_rate * 100
  )


K_Meh3_economic_df <- K_Meh3_economic_df %>%
  mutate(
    incorrect_acceptance_class = case_when(
      incorrect_acceptance_percent <= 2.5 ~ "0–2.5",
      incorrect_acceptance_percent <= 5   ~ "2.5–5",
      incorrect_acceptance_percent <= 7.5 ~ "5–7.5",
      incorrect_acceptance_percent <= 10  ~ "7.5–10",
      TRUE                                ~ ">10"
    )
  )

K_Meh3_economic_df <- K_Meh3_economic_df %>%
  mutate(
    incorrect_acceptance_class = factor(
      incorrect_acceptance_class,
      levels = c("0–2.5", "2.5–5", "5–7.5", "7.5–10", ">10")
    )
  )

K_Meh3_economic_df$title <- "bold(K~(mg~kg^{-1})*' - TabPFN')"

K_Meh3_Incorrect_Acceptance_Plot <- ggplot(
  K_Meh3_economic_df,
  aes(
    x = delta0,
    y = b,
    fill = incorrect_acceptance_class
  )
) +
  geom_raster() +

  scale_fill_manual(
    values = c(
      "0–2.5" = "#1F4E79",
      "2.5–5" = "#B2E2E2",
      "5–7.5" = "#FEC44F",
      "7.5–10" = "#FB6A4A",
      ">10" = "#67000D"
    ),
    name = expression("Incorrect Acceptances Among Accepted Predictions (%)" * ", " * alpha == 0.05),
    drop = FALSE
  ) +

  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(
    xlim = c(0,80),
    ylim = c(0,0.4),
    expand = FALSE
  )+

  labs(
    x =   expression(delta[0] ~ "(absolute threshold)"),
    y = expression(italic(b) ~ "(threshold scaling)")
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.title.y = element_text(size = 20, hjust = -0.7,vjust=1.9),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    legend.text = element_text(size = 18, colour = "black"),
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
    axis.line = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 18, colour = "black", face = "bold"),
    strip.background = element_rect(fill = "grey95", linewidth = 0.5),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 10)
  ) +
  facet_grid(. ~ title, labeller = label_parsed)

K_Meh3_Incorrect_Acceptance_Plot



#----------------------------------------------------------------------------------------------------------------
# P


P_Meh3_pred_values     <- as.numeric(P_Meh3_pred$P_Meh3_pred)
P_Meh3_true_values     <- as.numeric(P_Meh3_true$P_Meh3_target)
P_Meh3_quantile_matrix <- as.matrix(P_Meh3_pred_quantiles)

P_Meh3_error <- abs(P_Meh3_pred_values - P_Meh3_true_values)

delta0_values <- seq(-0.5, 60.5, length.out = 100)
b_values      <- seq(-0.01, 0.41, length.out = 100)

P_Meh3_economic_df <- expand.grid(
  delta0 = delta0_values,
  b      = b_values
) %>%
  arrange(b, delta0)

alpha_fixed <- 0.05




P_Meh3_economic_df$incorrect_acceptance_rate <- future_sapply(
  seq_len(nrow(P_Meh3_economic_df)),
  function(i) {

    a_i <- P_Meh3_economic_df$delta0[i]
    b_i <- P_Meh3_economic_df$b[i]

    p_vals <- compute_p_exceed_hybrid(
      quantile_matrix = P_Meh3_quantile_matrix,
      pred_values     = P_Meh3_pred_values,
      a               = a_i,
      b               = b_i,
      quantile_probs  = quantiles
    )

    delta_true <- a_i + b_i * P_Meh3_true_values

    accepted <- p_vals$p_value < alpha_fixed
    wrong    <- P_Meh3_error > delta_true

    if (sum(accepted) == 0) {
      0
    } else {
      mean(wrong[accepted])
    }

  },
  future.seed = TRUE
)

P_Meh3_economic_df <- P_Meh3_economic_df %>%
  mutate(
    incorrect_acceptance_percent = incorrect_acceptance_rate * 100
  )

P_Meh3_economic_df <- P_Meh3_economic_df %>%
  mutate(
    incorrect_acceptance_class = case_when(
      incorrect_acceptance_percent <= 2.5 ~ "0–2.5",
      incorrect_acceptance_percent <= 5   ~ "2.5–5",
      incorrect_acceptance_percent <= 7.5 ~ "5–7.5",
      incorrect_acceptance_percent <= 10  ~ "7.5–10",
      TRUE                                ~ ">10"
    )
  )

P_Meh3_economic_df <- P_Meh3_economic_df %>%
  mutate(
    incorrect_acceptance_class = factor(
      incorrect_acceptance_class,
      levels = c("0–2.5", "2.5–5", "5–7.5", "7.5–10", ">10")
    )
  )

P_Meh3_economic_df$title <- "bold(P~(mg~kg^{-1})*' - TabPFN')"

P_Meh3_Incorrect_Acceptance_Plot <- ggplot(
  P_Meh3_economic_df,
  aes(
    x = delta0,
    y = b,
    fill = incorrect_acceptance_class
  )
) +
  geom_raster() +
  scale_fill_manual(
    values = c(
      "0–2.5" = "#1F4E79",
      "2.5–5" = "#B2E2E2",
      "5–7.5" = "#FEC44F",
      "7.5–10" = "#FB6A4A",
      ">10" = "#67000D"
    ),
    name = expression("Incorrect Acceptances Among Accepted Predictions (%)" * ", " * alpha == 0.05),
    drop = FALSE
  ) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(
    xlim = c(0,60),
    ylim = c(0,0.4),
    expand = FALSE
  ) +
  labs(
    x = expression(delta[0]),
    y = expression(italic(b))
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 14, colour = "black"),
    legend.text = element_text(size = 18, colour = "black"),
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
    axis.line = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 18, colour = "black", face = "bold"),
    strip.background = element_rect(fill = "grey95", linewidth = 0.5)
  ) +
  facet_grid(. ~ title, labeller = label_parsed)

P_Meh3_Incorrect_Acceptance_Plot



# Rejector_Matrix_Incorrect_Acceptance_Plot_Map <- egg::ggarrange(
#   Clay_Incorrect_Acceptance_Plot,
#   SOM_Incorrect_Acceptance_Plot,
#   K_Meh3_Incorrect_Acceptance_Plot,
#   P_Meh3_Incorrect_Acceptance_Plot,
#   ncol=2,
#   nrow=2
# )
#
# Rejector_Matrix_Incorrect_Acceptance_Plot_Map
#
# ggsave(
#   "Figures/Rejector_Matrix_Incorrect_Acceptance_Plot_Map.jpg",
#   plot = Rejector_Matrix_Incorrect_Acceptance_Plot_Map,
#   width = 3000,
#   height = 2600,
#   units = "px",
#   dpi = 200
# )


# Additional Plots:


Total_C_true<- read.csv(file = "Results/TabPFN/Outer/Total_C_true.csv")
Sand_true<- read.csv(file = "Results/TabICL/Outer/Sand_true.csv")
Al_Meh3_true<-  read.csv(file = "Results/TabICL/Outer/Al_Meh3_true.csv")
pH_SMP_true<-  read.csv(file = "Results/TabPFN/Outer/pH_SMP_true.csv")


Total_C_pred<- read.csv(file = "Results/TabPFN/Outer/Total_C_pred.csv")
Sand_pred<- read.csv(file = "Results/TabICL/Outer/Sand_pred.csv")
Al_Meh3_pred<-  read.csv(file = "Results/TabICL/Outer/Al_Meh3_pred.csv")
pH_SMP_pred<-  read.csv(file = "Results/TabPFN/Outer/pH_SMP_pred.csv")


Total_C_pred_quantiles<- read.csv(file = "Results/TabPFN/Outer/Total_C_pred_quantiles.csv")
Sand_pred_quantiles<- read.csv(file = "Results/TabICL/Outer/Sand_pred_quantiles.csv")
Al_Meh3_pred_quantiles<-  read.csv(file = "Results/TabICL/Outer/Al_Meh3_pred_quantiles.csv")
pH_SMP_pred_quantiles<-  read.csv(file = "Results/TabPFN/Outer/pH_SMP_pred_quantiles.csv")


#----------------------------------------------------------------------------------------------------------------
# Sand



Sand_pred_values      <- as.numeric(Sand_pred$Sand_pred) * 10
Sand_true_values      <- as.numeric(Sand_true$Sand_target) * 10
Sand_quantile_matrix  <- as.matrix(Sand_pred_quantiles) * 10

Sand_error <- abs(Sand_pred_values - Sand_true_values)

delta0_values <- seq(-2, 172, length.out = 100)
b_values      <- seq(-0.11, 0.11, length.out = 100)
b_values

Sand_economic_df <- expand.grid(
  delta0 = delta0_values,
  b      = b_values
) %>%
  arrange(b, delta0)

Sand_economic_df$incorrect_acceptance_rate <- NA_real_

alpha_fixed <- 0.05

Sand_economic_df$incorrect_acceptance_rate <- future_sapply(
  seq_len(nrow(Sand_economic_df)),
  function(i) {

    a_i <- Sand_economic_df$delta0[i]
    b_i <- Sand_economic_df$b[i]

    p_vals <- compute_p_exceed_hybrid(
      quantile_matrix = Sand_quantile_matrix,
      pred_values     = Sand_pred_values,
      a               = a_i,
      b               = b_i,
      quantile_probs  = quantiles
    )

    delta_true <- a_i + b_i * Sand_true_values

    accepted <- p_vals$p_value < alpha_fixed
    wrong    <- Sand_error > delta_true

    if (sum(accepted) == 0) {
      0
    } else {
      mean(wrong[accepted])
    }

  },
  future.seed = TRUE
)

Sand_economic_df <- Sand_economic_df %>%
  mutate(
    incorrect_acceptance_percent = incorrect_acceptance_rate * 100
  )



Sand_economic_df <- Sand_economic_df %>%
  mutate(
    incorrect_acceptance_class = case_when(
      incorrect_acceptance_percent <= 2.5 ~ "0–2.5",
      incorrect_acceptance_percent <= 5   ~ "2.5–5",
      incorrect_acceptance_percent <= 7.5 ~ "5–7.5",
      incorrect_acceptance_percent <= 10  ~ "7.5–10",
      TRUE                                ~ ">10"
    )
  )

Sand_economic_df <- Sand_economic_df %>%
  mutate(
    incorrect_acceptance_class = factor(
      incorrect_acceptance_class,
      levels = c("0–2.5", "2.5–5", "5–7.5", "7.5–10", ">10")
    )
  )


Sand_economic_df$title <- "bold(Sand~(g~kg^{-1})*' - TabICL')"


Sand_Incorrect_Acceptance_Plot <- ggplot(
  Sand_economic_df,
  aes(
    x = delta0,
    y = b,
    fill = incorrect_acceptance_class
  )
) +

  geom_raster() +

  scale_fill_manual(
    values = c(
      "0–2.5" = "#1F4E79",
      "2.5–5" = "#B2E2E2",
      "5–7.5" = "#FEC44F",
      "7.5–10" = "#FB6A4A",
      ">10" = "#67000D"
    ),
    name = expression("Incorrect Acceptances Among Accepted Predictions (%)" * ", " * alpha == 0.05),
    drop = FALSE
  ) +

  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(
    xlim = c(0,170),
    ylim = c(-0.1,0.1),
    expand = FALSE
  )+

  labs(
    x = expression(delta[0]),
    y = expression(italic(b))
  ) +

  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.title = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    legend.text = element_text(size = 18, colour = "black"),
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
    axis.line = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 18, colour = "black", face = "bold"),
    strip.background = element_rect(fill = "grey95", linewidth = 0.5)
  ) +
  facet_grid(. ~ title, labeller = label_parsed)

Sand_Incorrect_Acceptance_Plot




#----------------------------------------------------------------------------------------------------------------
# Total_C



Total_C_pred_values      <- as.numeric(Total_C_pred$Total_C_pred) * 10
Total_C_true_values      <- as.numeric(Total_C_true$Total_C_target) * 10
Total_C_quantile_matrix  <- as.matrix(Total_C_pred_quantiles) * 10

Total_C_error <- abs(Total_C_pred_values - Total_C_true_values)

delta0_values <- seq(-0.05, 25.5, length.out = 100)
b_values      <- seq(-0.11, 0.11, length.out = 100)


Total_C_economic_df <- expand.grid(
  delta0 = delta0_values,
  b      = b_values
) %>%
  arrange(b, delta0)

Total_C_economic_df$incorrect_acceptance_rate <- NA_real_

alpha_fixed <- 0.05

Total_C_economic_df$incorrect_acceptance_rate <- future_sapply(
  seq_len(nrow(Total_C_economic_df)),
  function(i) {

    a_i <- Total_C_economic_df$delta0[i]
    b_i <- Total_C_economic_df$b[i]

    p_vals <- compute_p_exceed_hybrid(
      quantile_matrix = Total_C_quantile_matrix,
      pred_values     = Total_C_pred_values,
      a               = a_i,
      b               = b_i,
      quantile_probs  = quantiles
    )

    delta_true <- a_i + b_i * Total_C_true_values

    accepted <- p_vals$p_value < alpha_fixed
    wrong    <- Total_C_error > delta_true

    if (sum(accepted) == 0) {
      0
    } else {
      mean(wrong[accepted])
    }

  },
  future.seed = TRUE
)

Total_C_economic_df <- Total_C_economic_df %>%
  mutate(
    incorrect_acceptance_percent = incorrect_acceptance_rate * 100
  )



Total_C_economic_df <- Total_C_economic_df %>%
  mutate(
    incorrect_acceptance_class = case_when(
      incorrect_acceptance_percent <= 2.5 ~ "0–2.5",
      incorrect_acceptance_percent <= 5   ~ "2.5–5",
      incorrect_acceptance_percent <= 7.5 ~ "5–7.5",
      incorrect_acceptance_percent <= 10  ~ "7.5–10",
      TRUE                                ~ ">10"
    )
  )

Total_C_economic_df <- Total_C_economic_df %>%
  mutate(
    incorrect_acceptance_class = factor(
      incorrect_acceptance_class,
      levels = c("0–2.5", "2.5–5", "5–7.5", "7.5–10", ">10")
    )
  )


Total_C_economic_df$title <- "bold(TC~(g~kg^{-1})*' - TabPFN')"


Total_C_Incorrect_Acceptance_Plot <- ggplot(
  Total_C_economic_df,
  aes(
    x = delta0,
    y = b,
    fill = incorrect_acceptance_class
  )
) +

  geom_raster() +

  scale_fill_manual(
    values = c(
      "0–2.5" = "#1F4E79",
      "2.5–5" = "#B2E2E2",
      "5–7.5" = "#FEC44F",
      "7.5–10" = "#FB6A4A",
      ">10" = "#67000D"
    ),
    name = expression("Incorrect Acceptances Among Accepted Predictions (%)" * ", " * alpha == 0.05),
    drop = FALSE
  ) +

  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(
    xlim = c(0,25),
    ylim = c(-0.1,0.1),
    expand = FALSE
  )+

  labs(
    x = expression(delta[0]),
    y = expression(italic(b))
  ) +

  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 14, colour = "black"),
    legend.text = element_text(size = 18, colour = "black"),
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
    axis.line = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 18, colour = "black", face = "bold"),
    strip.background = element_rect(fill = "grey95", linewidth = 0.5)
  ) +
  facet_grid(. ~ title, labeller = label_parsed)

Total_C_Incorrect_Acceptance_Plot




#----------------------------------------------------------------------------------------------------------------
# Al





Al_Meh3_pred_values     <- as.numeric(Al_Meh3_pred$Al_Meh3_pred)
Al_Meh3_true_values     <- as.numeric(Al_Meh3_true$Al_Meh3_target)
Al_Meh3_quantile_matrix <- as.matrix(Al_Meh3_pred_quantiles)

Al_Meh3_error <- abs(Al_Meh3_pred_values - Al_Meh3_true_values)

delta0_values <- seq(-1, 501, length.out = 100)
b_values      <- seq(-0.11, 0.11, length.out = 100)

Al_Meh3_economic_df <- expand.grid(
  delta0 = delta0_values,
  b      = b_values
) %>%
  arrange(b, delta0)



alpha_fixed <- 0.05

Al_Meh3_economic_df$incorrect_acceptance_rate <- future_sapply(
  seq_len(nrow(Al_Meh3_economic_df)),
  function(i) {

    a_i <- Al_Meh3_economic_df$delta0[i]
    b_i <- Al_Meh3_economic_df$b[i]

    p_vals <- compute_p_exceed_hybrid(
      quantile_matrix = Al_Meh3_quantile_matrix,
      pred_values     = Al_Meh3_pred_values,
      a               = a_i,
      b               = b_i,
      quantile_probs  = quantiles
    )

    delta_true <- a_i + b_i * Al_Meh3_true_values

    accepted <- p_vals$p_value < alpha_fixed
    wrong    <- Al_Meh3_error > delta_true

    if (sum(accepted) == 0) {
      0
    } else {
      mean(wrong[accepted])
    }

  },
  future.seed = TRUE
)

Al_Meh3_economic_df <- Al_Meh3_economic_df %>%
  mutate(
    incorrect_acceptance_percent = incorrect_acceptance_rate * 100
  )


Al_Meh3_economic_df <- Al_Meh3_economic_df %>%
  mutate(
    incorrect_acceptance_class = case_when(
      incorrect_acceptance_percent <= 2.5 ~ "0–2.5",
      incorrect_acceptance_percent <= 5   ~ "2.5–5",
      incorrect_acceptance_percent <= 7.5 ~ "5–7.5",
      incorrect_acceptance_percent <= 10  ~ "7.5–10",
      TRUE                                ~ ">10"
    )
  )

Al_Meh3_economic_df <- Al_Meh3_economic_df %>%
  mutate(
    incorrect_acceptance_class = factor(
      incorrect_acceptance_class,
      levels = c("0–2.5", "2.5–5", "5–7.5", "7.5–10", ">10")
    )
  )

Al_Meh3_economic_df$title <- "bold(Al~(mg~kg^{-1})*' - TabICL')"


Al_Meh3_Incorrect_Acceptance_Plot <- ggplot(
  Al_Meh3_economic_df,
  aes(
    x = delta0,
    y = b,
    fill = incorrect_acceptance_class
  )
) +
  geom_raster() +

  scale_fill_manual(
    values = c(
      "0–2.5" = "#1F4E79",
      "2.5–5" = "#B2E2E2",
      "5–7.5" = "#FEC44F",
      "7.5–10" = "#FB6A4A",
      ">10" = "#67000D"
    ),
    name = expression("Incorrect Acceptances Among Accepted Predictions (%)" * ", " * alpha == 0.05),
    drop = FALSE
  ) +

  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(
    xlim = c(0,500),
    ylim = c(-0.1,0.1),
    expand = FALSE
  )+

  labs(
    x =   expression(delta[0] ~ "(absolute threshold)"),
    y = expression(italic(b))
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.title.x = element_text(size = 20,hjust = 1.3, face = "italic"),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    legend.text = element_text(size = 18, colour = "black"),
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
    axis.line = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 18, colour = "black", face = "bold"),
    strip.background = element_rect(fill = "grey95", linewidth = 0.5)
  ) +
  facet_grid(. ~ title, labeller = label_parsed)

Al_Meh3_Incorrect_Acceptance_Plot


#----------------------------------------------------------------------------------------------------------------
# pH



pH_SMP_pred_values     <- as.numeric(pH_SMP_pred$pH_SMP_pred)
pH_SMP_true_values     <- as.numeric(pH_SMP_true$pH_SMP_target)
pH_SMP_quantile_matrix <- as.matrix(pH_SMP_pred_quantiles)

pH_SMP_error <- abs(pH_SMP_pred_values - pH_SMP_true_values)

delta0_values <- seq(0.01, 0.5, length.out = 100)
b_values      <- seq(-0.11, 0.11, length.out = 100)

pH_SMP_economic_df <- expand.grid(
  delta0 = delta0_values,
  b      = b_values
) %>%
  arrange(b, delta0)

alpH_SMPa_fixed <- 0.05

pH_SMP_economic_df

pH_SMP_economic_df$incorrect_acceptance_rate <- future_sapply(
  seq_len(nrow(pH_SMP_economic_df)),
  function(i) {

    a_i <- pH_SMP_economic_df$delta0[i]
    b_i <- pH_SMP_economic_df$b[i]

    p_vals <- compute_p_exceed_hybrid(
      quantile_matrix = pH_SMP_quantile_matrix,
      pred_values     = pH_SMP_pred_values,
      a               = a_i,
      b               = b_i,
      quantile_probs  = quantiles
    )

    delta_true <- a_i + b_i * pH_SMP_true_values

    accepted <- p_vals$p_value < alpH_SMPa_fixed
    wrong    <- pH_SMP_error > delta_true

    if (sum(accepted) == 0) {
      0
    } else {
      mean(wrong[accepted])
    }

  },
  future.seed = TRUE
)


pH_SMP_economic_df <- pH_SMP_economic_df %>%
  mutate(
    incorrect_acceptance_percent = incorrect_acceptance_rate * 100
  )

pH_SMP_economic_df <- pH_SMP_economic_df %>%
  mutate(
    incorrect_acceptance_class = case_when(
      incorrect_acceptance_percent <= 2.5 ~ "0–2.5",
      incorrect_acceptance_percent <= 5   ~ "2.5–5",
      incorrect_acceptance_percent <= 7.5 ~ "5–7.5",
      incorrect_acceptance_percent <= 10  ~ "7.5–10",
      TRUE                                ~ ">10"
    )
  )

pH_SMP_economic_df <- pH_SMP_economic_df %>%
  mutate(
    incorrect_acceptance_class = factor(
      incorrect_acceptance_class,
      levels = c("0–2.5", "2.5–5", "5–7.5", "7.5–10", ">10")
    )
  )

pH_SMP_economic_df$title <- "bold(pH*' - TabPFN')"

pH_SMP_Incorrect_Acceptance_Plot <- ggplot(
  pH_SMP_economic_df,
  aes(
    x = delta0,
    y = b,
    fill = incorrect_acceptance_class
  )
) +
  geom_raster() +
  scale_fill_manual(
    values = c(
      "0–2.5" = "#1F4E79",
      "2.5–5" = "#B2E2E2",
      "5–7.5" = "#FEC44F",
      "7.5–10" = "#FB6A4A",
      ">10" = "#67000D"
    ),
    name = expression("Incorrect Acceptances Among Accepted Predictions (%)" * ", " * alpH_SMPa == 0.05),
    drop = FALSE
  ) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(
    xlim = c(0.025, 0.5),
    ylim = c(-0.1,0.1),
    expand = FALSE
  )+
  labs(
    x = expression(delta[0]),
    y = expression(italic(b))
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 14, colour = "black"),
    legend.text = element_text(size = 18, colour = "black"),
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
    axis.line = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 18, colour = "black", face = "bold"),
    strip.background = element_rect(fill = "grey95", linewidth = 0.5)
  ) +
  facet_grid(. ~ title, labeller = label_parsed)

pH_SMP_Incorrect_Acceptance_Plot










Rejector_Matrix_Incorrect_Acceptance_Plot_Map_8 <- egg::ggarrange(
  Clay_Incorrect_Acceptance_Plot,
  SOM_Incorrect_Acceptance_Plot,
  K_Meh3_Incorrect_Acceptance_Plot,
  P_Meh3_Incorrect_Acceptance_Plot,
  Sand_Incorrect_Acceptance_Plot,
  Total_C_Incorrect_Acceptance_Plot,
  Al_Meh3_Incorrect_Acceptance_Plot,
  pH_SMP_Incorrect_Acceptance_Plot,
  ncol=2,
  nrow=4
)


ggsave(
  "Figures/FigS5.jpg",
  plot = Rejector_Matrix_Incorrect_Acceptance_Plot_Map_8,
  width = 3000*1.5,
  height = 4500*1.5,
  units = "px",
  dpi = 300
)
