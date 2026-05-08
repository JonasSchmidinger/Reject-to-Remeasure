# In this script, we create Fig.S7

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

plan(multisession, workers = 6)



quantiles <- c(
  0.001,
  round(seq(0.005, 0.995, by = 0.005), 3),
  0.999
)


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
      x = c(q_vals),
      y = c(0, q_probs),
      right = TRUE
    )

    y_hat <- pred_values[i]

    #Need different criteria not based on y_hat but on Y, as changes along Y will also change the threshold (as compared to absolute threshold with fixed threshold widths)
    L <- (y_hat - a) / (1 + b)
    U <- (y_hat + a) / (1 - b)

    p_vals[i] <- cdf_fun(L) + (1 - cdf_fun(U))
  }

  data.frame(p_value = p_vals)
}

compute_p_exceed_absolute <- function(quantile_matrix,
                                      pred_values,
                                      a,
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

    # Absolute bounds
    L <- y_hat - a
    U <- y_hat + a

    # Probability of bad event
    p_vals[i] <- cdf_fun(L) + (1 - cdf_fun(U))
  }

  data.frame(p_value = p_vals)
}




#----------------------------------------------------------------------------------------------------------------
# Sand

Sand_pred_values      <- as.numeric(Sand_pred$Sand_pred) * 10
Sand_quantile_matrix  <- as.matrix(Sand_pred_quantiles) * 10

delta0_values <- seq(-2, 172, by = 2)
b_values      <- seq(-0.11, 0.11, by = 0.005)
b_values


Sand_economic_df <- expand.grid(
  delta0 = delta0_values,
  b      = b_values
) %>%
  arrange(b, delta0)

Sand_economic_df$acceptance_rate <- NA_real_

alpha_fixed <- 0.05


Sand_economic_df$acceptance_rate <- future_sapply(
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

    mean(p_vals$p_value < alpha_fixed)

  },
  future.seed = TRUE
)

Sand_economic_df <- Sand_economic_df %>%
  mutate(
    acceptance_percent = acceptance_rate * 100
  )


Sand_economic_df$title <- "bold(Sand~(g~kg^{-1})*' - TabICL')"

contour_levels <- seq(10, 90, by = 20)
contour_levels

Sand_economic_df

Sand_Economic_Plot <- ggplot(
  Sand_economic_df,
  aes(
    x = delta0,
    y = b,
    fill = acceptance_percent
  )
) +

  geom_raster() +

  geom_contour(
    aes(z = acceptance_percent),
    breaks = contour_levels,
    colour = "black",
    linewidth = 0.4
  ) +

  geom_text_contour(
    aes(z = acceptance_percent),
    breaks = contour_levels,
    size = 6,
    stroke = 0.1,
    family = "Arial",
    skip = 0,
    rotate = TRUE
  )+

  scale_fill_gradientn(
    colours = c(
      "#67000D",
      "darkred",
      "#FB6A4A",
      "#FEC44F",
      "#B2E2E2",
      "#6BAED6",
      "#1F4E79"
    ),
    values = scales::rescale(
      c(1, 3, 10, 30, 55, 80, 100)
    ),
    limits = c(0, 100),
    name = "Acceptance (%)"
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
    y = expression(italic(b) ~ "(threshold scaling)")
  ) +

  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.title.y = element_text(size = 20, hjust = -0.8,vjust= 3),
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

Sand_Economic_Plot



#----------------------------------------------------------------------------------------------------------------
# Total_C


Total_C_pred_values     <- as.numeric(Total_C_pred$Total_C_pred) * 10
Total_C_quantile_matrix <- as.matrix(Total_C_pred_quantiles) * 10

delta0_values <- seq(-0.05, 25.5, by = 0.5)
b_values      <- seq(-0.11, 0.11, by = 0.005)


Total_C_economic_df <- expand.grid(
  delta0 = delta0_values,
  b      = b_values
) %>%
  arrange(b, delta0)

alpha_fixed <- 0.05

Total_C_economic_df$acceptance_rate <- future_sapply(
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

    mean(p_vals$p_value < alpha_fixed)

  },
  future.seed = TRUE
)

Total_C_economic_df <- Total_C_economic_df %>%
  mutate(
    acceptance_percent = acceptance_rate * 100
  )

Total_C_economic_df$title <- "bold(TC~(g~kg^{-1})*' - TabPFN')"

Total_C_Economic_Plot <- ggplot(
  Total_C_economic_df,
  aes(
    x = delta0,
    y = b,
    fill = acceptance_percent
  )
) +
  geom_raster() +

  geom_contour(
    aes(z = acceptance_percent),
    breaks = contour_levels,
    colour = "black",
    linewidth = 0.4
  ) +

  geom_text_contour(
    aes(z = acceptance_percent),
    breaks = contour_levels,
    size = 6,
    stroke = 0.1,
    family = "Arial",
    skip = 0,
    rotate = TRUE
  )+

  scale_fill_gradientn(
    colours = c(
      "#67000D",
      "darkred",
      "#FB6A4A",
      "#FEC44F",
      "#B2E2E2",
      "#6BAED6",
      "#1F4E79"
    ),
    values = scales::rescale(c(1, 3, 10, 30, 55, 80, 100)),
    limits = c(0, 100),
    name = expression("Acceptance Rate (%)" * ", " * alpha == 0.05)  ) +

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
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    legend.position = c(-0.35, 1.1),
    legend.justification = c(0, 0),
    plot.margin = margin(t = 85, r = 10, b = 10, l = 10),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(size = 20, colour = "black"),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.text = element_text(size = 18, colour = "black", margin = margin(t = 3, r = 25)),
    legend.ticks = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
    axis.line = element_blank(),
    strip.text = element_text(size = 18, colour = "black", face = "bold"),
    strip.background = element_rect(fill = "grey95", linewidth = 0.5)
  ) +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(12, "cm"),
      barheight = unit(1, "cm"),
      ticks = TRUE,
      frame.colour = "black",
      order = 1
    )
  ) +
  facet_grid(. ~ title, labeller = label_parsed)

Total_C_Economic_Plot



#----------------------------------------------------------------------------------------------------------------
# Al


Al_Meh3_pred_values     <- as.numeric(Al_Meh3_pred$Al_Meh3_pred)
Al_Meh3_quantile_matrix <- as.matrix(Al_Meh3_pred_quantiles)

delta0_values <- seq(-1, 501, by = 5)
b_values      <- seq(-0.11, 0.11, by = 0.005)

Al_Meh3_economic_df <- expand.grid(
  delta0 = delta0_values,
  b      = b_values
) %>%
  arrange(b, delta0)

alpha_fixed <- 0.05

Al_Meh3_economic_df$acceptance_rate <- future_sapply(
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

    mean(p_vals$p_value < alpha_fixed)

  },
  future.seed = TRUE
)

Al_Meh3_economic_df <- Al_Meh3_economic_df %>%
  mutate(
    acceptance_percent = acceptance_rate * 100
  )

Al_Meh3_economic_df$title <- "bold(Al~(mg~kg^{-1})*' - TabICL')"

Al_Meh3_Economic_Plot <- ggplot(
  Al_Meh3_economic_df,
  aes(
    x = delta0,
    y = b,
    fill = acceptance_percent
  )
) +
  geom_raster() +

  geom_contour(
    aes(z = acceptance_percent),
    breaks = contour_levels,
    colour = "black",
    linewidth = 0.4
  ) +

  geom_text_contour(
    aes(z = acceptance_percent),
    breaks = contour_levels,
    size = 6,
    stroke = 0.1,
    family = "Arial",
    skip = 0,
    rotate = TRUE
  )+

  scale_fill_gradientn(
    colours = c(
      "#67000D",
      "darkred",
      "#FB6A4A",
      "#FEC44F",
      "#B2E2E2",
      "#6BAED6",
      "#1F4E79"
    ),
    values = scales::rescale(c(1,3,10,30,55,80,100)),
    limits = c(0,100)
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
    axis.title.x = element_text(size = 20,hjust = 1.35, face = "italic"),
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
Al_Meh3_Economic_Plot



#----------------------------------------------------------------------------------------------------------------
# pH




pH_SMP_pred_values     <- as.numeric(pH_SMP_pred$pH_SMP_pred)
pH_SMP_quantile_matrix <- as.matrix(pH_SMP_pred_quantiles)

delta0_values <- seq(0, 0.5, by = 0.005)
b_values      <- seq(-0.11, 0.11, by = 0.005)

pH_SMP_economic_df <- expand.grid(
  delta0 = delta0_values,
  b      = b_values
) %>%
  arrange(b, delta0)

alpha_fixed <- 0.05

pH_SMP_economic_df$acceptance_rate <- future_sapply(
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

    mean(p_vals$p_value < alpha_fixed)

  },
  future.seed = TRUE
)

pH_SMP_economic_df <- pH_SMP_economic_df %>%
  mutate(
    acceptance_percent = acceptance_rate * 100
  )

pH_SMP_economic_df$title <-"bold(pH*' - TabPFN')"

pH_SMP_Economic_Plot <- ggplot(
  pH_SMP_economic_df,
  aes(
    x = delta0,
    y = b,
    fill = acceptance_percent
  )
) +
  geom_raster() +

  geom_contour(
    aes(z = acceptance_percent),
    breaks = contour_levels,
    colour = "black",
    linewidth = 0.4
  ) +

  geom_text_contour(
    aes(z = acceptance_percent),
    breaks = contour_levels,
    size = 6,
    stroke = 0.1,
    family = "Arial",
    skip = 0,
    rotate = TRUE
  )+

  scale_fill_gradientn(
    colours = c(
      "#67000D",
      "darkred",
      "#FB6A4A",
      "#FEC44F",
      "#B2E2E2",
      "#6BAED6",
      "#1F4E79"
    ),
    values = scales::rescale(c(1,3,10,30,55,80,100)),
    limits = c(0,100)
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
    axis.text.y =  element_blank(),
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

pH_SMP_Economic_Plot





Rejector_Matrix_Economic_Plot_Map_Appendix <- egg::ggarrange(
  Sand_Economic_Plot,
  Total_C_Economic_Plot,
  Al_Meh3_Economic_Plot,
  pH_SMP_Economic_Plot, ncol=2,nrow=2)

Rejector_Matrix_Economic_Plot_Map_Appendix

ggsave("Figures/FigS7.jpg", plot = Rejector_Matrix_Economic_Plot_Map_Appendix, width = 3000 * 1.5, height = 2600 * 1.5, units = "px", dpi = 300)


