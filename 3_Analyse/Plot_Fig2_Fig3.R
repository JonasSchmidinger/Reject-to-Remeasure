# In this script, we create Fig.2 and Fig.3

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
library(showtext)
library(sysfonts)

font_add(
  family = "Arial",
  regular = "C:/Windows/Fonts/arial.ttf",
  italic  = "C:/Windows/Fonts/ariali.ttf",
  bold    = "C:/Windows/Fonts/arialbd.ttf",
  bolditalic = "C:/Windows/Fonts/arialbi.ttf"
)


sym_yhat  <- "\u0177"  # ŷ
sym_delta <- "\u03B4"  # δ
sym_minus <- "\u2212"  # − (proper minus)
sym_geq   <- "\u2265"  # ≥
sym_mid   <- "\u2223"



quantiles <-
  c(
    0.001,
    round(seq(0.005, 0.995, by = 0.005), 3),
    0.999
  )
quantiles


calculate_mae_rejection <- function(error, is_rejected) {
  c(
    MAE_before_rejection = mean(error, na.rm = TRUE),
    MAE_after_rejection  = mean(error[!is_rejected], na.rm = TRUE)
  )
}


SOM_true<- read.csv(file = "Results/TabPFN/Outer/SOM_true.csv")
Clay_true<- read.csv(file = "Results/TabPFN/Outer/Clay_true.csv")
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





compute_p_exceed_all <- function(quantile_matrix, pred_values, threshold, quantile_probs) {

  n <- nrow(quantile_matrix)
  p_vals <- numeric(n)

  for (i in seq_len(n)) {

    q_vals <- as.numeric(quantile_matrix[i, ])

    ord <- order(q_vals)
    q_vals <- q_vals[ord]
    q_probs <- quantile_probs[ord]

    cdf_fun <- stepfun(
      x = c(q_vals),
      y = c(0, q_probs),
      right = TRUE
    )

    y_true <- as.numeric(pred_values[i, 1])

    p_vals[i] <- cdf_fun(y_true - threshold) +
      (1 - cdf_fun(y_true + threshold))
  }

  data.frame(p_value = p_vals)
}





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

    y_hat <- as.numeric(pred_values[i, 1])

    L <- (y_hat - a) / (1 + b)
    U <- (y_hat + a) / (1 - b)

    p_vals[i] <- cdf_fun(L) + (1 - cdf_fun(U))
  }

  data.frame(p_value = p_vals)
}




#----------------------------------------------------------------------------------------------------------------
# Clay



Clay_abs_buffer  <- 50
Clay_rel_factor  <- 0.15
Clay_threshold_alpha <- 0.05




Clay_p_value <- compute_p_exceed_hybrid(
  quantile_matrix = Clay_pred_quantiles * 10, # *10 because of changing from % to g kg^-1
  pred_values     = data.frame(pred = Clay_pred$Clay_pred * 10),
  a               = Clay_abs_buffer,
  b               = Clay_rel_factor,
  quantile_probs  = quantiles
)







Clay_error_p.df <- data.frame(
  pred    =  as.numeric(Clay_pred$Clay_pred) * 10,
  true    = as.numeric(Clay_true$Clay_target) * 10,
  p_value = Clay_p_value$p_value
) %>%
  mutate(
    pred  = as.numeric(pred),
    true  = as.numeric(true),
    error = abs(pred - true),

    # thresholds
    delta_model = Clay_abs_buffer + Clay_rel_factor * pred,  # for model-assumed probability
    delta_true  = Clay_abs_buffer + Clay_rel_factor * true,  # Real boundary for outliers

    # wrong means we exceeded the error threshold, rejected gives the actual guess of the rejector
    is_wrong    = error > delta_true,
    is_rejected = p_value > Clay_threshold_alpha
  )

Clay_error_p.df



Clay_confusion <- factor(
  ifelse(Clay_error_p.df$is_rejected &  Clay_error_p.df$is_wrong, "Correct Rejection",
         ifelse(Clay_error_p.df$is_rejected & !Clay_error_p.df$is_wrong, "Incorrect Rejection",
                ifelse(!Clay_error_p.df$is_rejected &  Clay_error_p.df$is_wrong, "Incorrect Acceptance",
                       "Correct Acceptance"))),
  levels = c(
    "Incorrect Rejection",
    "Incorrect Acceptance",
    "Correct Rejection",
    "Correct Acceptance"
  )
)

Clay_error_cor_plot.df <- Clay_error_p.df %>%
  mutate(
    confusion = Clay_confusion,
    title =  "bold(Clay~(g~kg^{-1})*' - TabICL')"
  )

confusion_labels_Clay <- Clay_error_cor_plot.df %>%
  mutate(
    group = ifelse(is_rejected, "Rejected", "Accepted")
  ) %>%
  count(confusion, group) %>%
  group_by(group) %>%
  mutate(
    pct_within = n / sum(n) * 100
  ) %>%
  ungroup() %>%
  mutate(
    pct_total = n / sum(n) * 100
  ) %>%
  mutate(
    label = sprintf(
      "%s: %.1f%% (%.1f%% of %s)",
      confusion,
      pct_total,
      pct_within,
      tolower(group)
    )
  )
confusion_labels_Clay

acceptance_rate <- Clay_error_cor_plot.df %>%
  summarise(rate = mean(!is_rejected) * 100) %>%
  pull(rate)
acceptance_rate

rejection_rate <- Clay_error_cor_plot.df %>%
  summarise(rate = mean(is_rejected) * 100) %>%
  pull(rate)
rejection_rate


Clay_error_cor_plot.df$p_value

Clay_Rejector_AE_alpha <-ggplot(
  Clay_error_cor_plot.df,
  aes(
    x = p_value,
    y = error,
    colour = confusion
  )
) +
  geom_point(
    shape = 21,
    aes(fill = confusion),
    colour = "black",
    stroke = 0.4,
    alpha = 0.5,
    size = 2.5
  )+
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  scale_fill_manual(
    values = c(
      "Incorrect Rejection" = "#67000D",
      "Incorrect Acceptance" = "#F4A6A6",
      "Correct Rejection"  = "#6BAED6",
      "Correct Acceptance"  = "#1F4E79"
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    y = "Absolute Error"
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.title.y = element_text(size = 20,hjust = -0.34),
    axis.text = element_text(size = 14,color="black"),
    axis.title.x =  element_blank(),
    legend.text = element_text(size = 18,color="black"),
    legend.title = element_blank(),
    strip.text = element_text(size = 18,colour = "black", face = "bold"),
    strip.background = element_rect(fill = "grey95",linewidth = 0.5 ),
    panel.border = element_rect(colour = "black", linewidth = 0.5,fill = NA ),
    axis.line = element_blank(),
    legend.position = "none"
  )+
  annotate(
    "text",
    x = 0.053,
    y = Inf,
    label = bquote(alpha == .(0.05)),
    vjust = 1.5,
    hjust = 0,
    size = 5,
    family = "Arial"
  )+

  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = confusion_labels_Clay$label,
    hjust = 1.02,
    vjust = seq(1.2, by = 1.3, length.out = nrow(confusion_labels_Clay)),
    colour = c(
      "Incorrect Rejection" = "#67000D",
      "Incorrect Acceptance" = "#F4A6A6",
      "Correct Rejection"  = "#6BAED6",
      "Correct Acceptance"  = "#1F4E79"
    )[confusion_labels_Clay$confusion],
    size = 5,
    family = "Arial"
  )+
  facet_grid(. ~ title, labeller = label_parsed)

Clay_Rejector_AE_alpha



Clay_MAE <- Clay_error_cor_plot.df %>%
  summarise(
    MAE_accepted = mean(error[!is_rejected], na.rm = TRUE),
    MAE_rejected = mean(error[is_rejected], na.rm = TRUE)
  )



Clay_error_cor_plot.df$title <- "bold(Clay~(g~kg^{-1})*' - TabICL')"

Clay_Rejector_1to1 <- ggplot(
  Clay_error_cor_plot.df,
  aes(x = true, y = pred, fill = p_value)
) +
  geom_point(
    shape = 21,
    colour = "black",
    stroke = 0.4,
    size = 2.5,
    alpha = 0.9
  ) +

  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_function(
    fun = function(x) x + (Clay_abs_buffer + Clay_rel_factor * x),
    linetype = "dotted",
    linewidth = 0.9
  ) +
  geom_function(
    fun = function(x) x - (Clay_abs_buffer + Clay_rel_factor * x),
    linetype = "dotted",
    linewidth = 0.9
  )+
  scale_fill_gradientn(
    colours = c(
      "#1F4E79","#6BAED6","#B2E2E2",
      "#FEC44F","#FB6A4A","darkred","#67000D"
    ),
    values = scales::rescale(
      c(0,0.015,0.05,0.0501,0.2,0.4,1),
      from = c(0,1)
    ),
    limits = c(0,1),
    breaks = c(0,0.05,1),
    labels = c("0","α = 0.05","1")
  ) +

  scale_x_continuous(
    limits = c(0, 1000),
    expand = c(0.005, 0.005)
  ) +
  scale_y_continuous(
    limits = c(0,  1000),
    expand = c(0.005, 0.005)
  ) +

  labs(
    x = "Measured Value",
    y = "Predicted Value"
  ) +

  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.title.y = element_text(size = 20,hjust = -0.6, vjust=2.3),
    axis.text = element_text(size = 14,color="black"),
    axis.title.x =  element_blank(),
    legend.text = element_text(size = 14,color="black"),
    legend.title = element_text(size = 18,color="black"),
    strip.text = element_text(size = 18,colour = "black", face = "bold"),
    strip.background = element_rect(fill = "grey95",linewidth = 0.5 ),
    panel.border = element_rect(colour = "black", linewidth = 0.5,fill = NA ),
    axis.line = element_blank(),
    legend.position = "none"
  )+
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = expression(
      italic("δ") * "(" * italic(Y) * ")" ==
        50 + 0.15 * italic(Y)
    ),
    parse = TRUE,
    hjust = -0.7,
    vjust = 1.7,
    size = 5,
    family = "Arial"
  )+
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = sprintf("MAE[rejected] == %.1f", Clay_MAE$MAE_rejected),
    parse = TRUE,
    hjust = 1.05,
    vjust = -2.2,
    colour = "#67000D",
    size = 5,
    family = "Arial"
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = sprintf("MAE[accepted] == %.1f", Clay_MAE$MAE_accepted),
    parse = TRUE,
    hjust = 1.05,
    vjust = -0.8,
    colour = "#1F4E79",
    size = 5,
    family = "Arial"
  ) +
  facet_grid(. ~ title, labeller = label_parsed)

Clay_Rejector_1to1



























#----------------------------------------------------------------------------------------------------------------
# SOM






SOM_abs_buffer  <- 2.5
SOM_rel_factor  <- 0.3
SOM_threshold_alpha <- 0.05




SOM_p_value <- compute_p_exceed_hybrid(
  quantile_matrix = SOM_pred_quantiles * 10,
  pred_values     = data.frame(pred = SOM_pred$SOM_pred * 10),
  a               = SOM_abs_buffer,
  b               = SOM_rel_factor,
  quantile_probs  = quantiles
)







SOM_error_p.df <- data.frame(
  pred    =  as.numeric(SOM_pred$SOM_pred) * 10,
  true    = as.numeric(SOM_true$SOM_target) * 10,
  p_value = SOM_p_value$p_value
) %>%
  mutate(
    pred  = as.numeric(pred),
    true  = as.numeric(true),
    error = abs(pred - true),

    delta_model = SOM_abs_buffer + SOM_rel_factor * pred,
    delta_true  = SOM_abs_buffer + SOM_rel_factor * true,

    is_wrong    = error > delta_true,
    is_rejected = p_value > SOM_threshold_alpha
  )

SOM_error_p.df



SOM_confusion <- factor(
  ifelse(SOM_error_p.df$is_rejected &  SOM_error_p.df$is_wrong, "Correct Rejection",
         ifelse(SOM_error_p.df$is_rejected & !SOM_error_p.df$is_wrong, "Incorrect Rejection",
                ifelse(!SOM_error_p.df$is_rejected &  SOM_error_p.df$is_wrong, "Incorrect Acceptance",
                       "Correct Acceptance"))),
  levels = c(
    "Incorrect Rejection",
    "Incorrect Acceptance",
    "Correct Rejection",
    "Correct Acceptance"
  )
)



SOM_error_cor_plot.df <- SOM_error_p.df %>%
  mutate(
    confusion = SOM_confusion,
    title = "bold(SOM~(g~kg^{-1})*' - TabPFN')"
  )

confusion_labels_SOM <- SOM_error_cor_plot.df %>%
  mutate(
    group = ifelse(is_rejected, "Rejected", "Accepted")
  ) %>%
  count(confusion, group) %>%
  group_by(group) %>%
  mutate(
    pct_within = n / sum(n) * 100
  ) %>%
  ungroup() %>%
  mutate(
    pct_total = n / sum(n) * 100
  ) %>%
  mutate(
    label = sprintf(
      "%s: %.1f%% (%.1f%% of %s)",
      confusion,
      pct_total,
      pct_within,
      tolower(group)
    )
  )
confusion_labels_SOM




SOM_error_cor_plot.df$p_value

SOM_Rejector_AE_alpha <-ggplot(
  SOM_error_cor_plot.df,
  aes(
    x = p_value,
    y = error ,
    colour = confusion
  )
) +
  geom_point(
    shape = 21,
    aes(fill = confusion),
    colour = "black",
    stroke = 0.4,
    alpha = 0.5,
    size = 2.5
  )+
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  scale_fill_manual(
    values = c(
      "Incorrect Rejection" = "#67000D",
      "Incorrect Acceptance" = "#F4A6A6",
      "Correct Rejection"  = "#6BAED6",
      "Correct Acceptance"  = "#1F4E79"
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    y = "Absolute Error"
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    legend.position = c(-1, 1.1),
    legend.justification = c(0, 0),
    plot.margin = margin(t = 85, r = 10, b = 10, l = 10),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.key.width = unit(1.5, "cm"),
    legend.ticks.length = unit(0.6, "cm"),
    legend.key.height = unit(1, "cm"),


    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    axis.text = element_text(size = 14,color="black"),
    axis.title =  element_blank(),
    legend.text = element_text(size = 18,color="black",margin = margin(r = 25)),
    strip.text = element_text(size = 18,colour = "black", face = "bold"),
    strip.background = element_rect(fill = "grey95",linewidth = 0.5 ),
    panel.border = element_rect(colour = "black", linewidth = 0.5,fill = NA ),
    axis.line = element_blank()
  )+
  guides(
    colour = "none",
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        size = 5,
        alpha = 0.8
      )
    )
  )+
  annotate(
    "text",
    x = 0.053,
    y = Inf,
    label = bquote(alpha == .(0.05)),
    vjust = 1.5,
    hjust = 0,
    size = 5,
    family = "Arial"
  )+
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = confusion_labels_SOM$label,
    hjust = 1.02,
    vjust = seq(1.2, by = 1.3, length.out = nrow(confusion_labels_SOM)),
    colour = c(
      "Incorrect Rejection" = "#67000D",
      "Incorrect Acceptance" = "#F4A6A6",
      "Correct Rejection"  = "#6BAED6",
      "Correct Acceptance"  = "#1F4E79"
    )[confusion_labels_SOM$confusion],
    size = 5,
    family = "Arial"
  )+
  facet_grid(. ~ title, labeller = label_parsed)

SOM_Rejector_AE_alpha


SOM_MAE <- SOM_error_cor_plot.df %>%
  summarise(
    MAE_accepted = mean(error[!is_rejected], na.rm = TRUE),
    MAE_rejected = mean(error[is_rejected], na.rm = TRUE)
  )



SOM_error_cor_plot.df$title <- "bold(SOM~(g~kg^{-1})*' - TabPFN')"
SOM_error_cor_plot.df

SOM_Rejector_1to1 <- ggplot(
  SOM_error_cor_plot.df,
  aes(x = true, y = pred, fill = p_value)) +
  geom_point(
    shape = 21,
    colour = "black",
    stroke = 0.4,
    size = 2.5,
    alpha = 0.9
  )+

  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 1
  ) +

  geom_function(
    fun = function(x) x + (SOM_abs_buffer + SOM_rel_factor * x),
    linetype = "dotted",
    linewidth = 0.9
  ) +
  geom_function(
    fun = function(x) x - (SOM_abs_buffer + SOM_rel_factor * x),
    linetype = "dotted",
    linewidth = 0.9
  )+

  scale_fill_gradientn(
    colours = c(
      "#1F4E79",
      "#6BAED6",
      "#B2E2E2",
      "#FEC44F",
      "#FB6A4A",
      "darkred",
      "#67000D"
    ),

    values = scales::rescale(
      c(
        0,
        0.015,
        0.05,
        0.0501,
        0.2,
        0.4,
        1
      ),
      from = c(0,1)
    ),

    limits = c(0,1),

    breaks = c(0,0.05,1),
    labels = c("0", expression(alpha == 0.05), "1")
  ) +
  guides(
    fill = guide_colourbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      label.hjust = 0,

      ticks = TRUE,
      ticks.colour = "black",
      ticks.linewidth = 0.7,

      barwidth  = unit(7.5, "cm"),
      barheight = unit(0.5, "cm"),

      frame.colour = "black",
      frame.linewidth = 0.9
    )
  )+
  scale_x_continuous(limits = c(0, 200),expand = c(0.005, 0.005)) +
  scale_y_continuous(limits = c(0, 200),expand = c(0.005, 0.005)) +

  labs(
    x = "Measured Value",
    y = "Predicted Value",
    fill = expression(
      plain("Pr") *
        "(" *
        "|" * italic(Y) * " − " * italic("ŷ") * "|" *
        italic(" ≥ ") *
        italic("δ") * "(" * italic(Y) * ")" *
        "|" * italic(x) *
        ")"
    )
  )+
  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    legend.position = c(0, 1.2),
    plot.margin = margin(t = 85, r = 10, b = 10, l = 10),
    legend.ticks.length = unit(0.5, "cm"),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    axis.text = element_text(size = 14,color="black"),
    axis.title =  element_blank(),
    legend.text = element_text(size = 14,color="black"),
    legend.title = element_text(size = 18,color="black"),
    strip.text = element_text(size = 18,colour = "black", face = "bold"),
    strip.background = element_rect(fill = "grey95",linewidth = 0.5 ),
    panel.border = element_rect(colour = "black", linewidth = 0.5,fill = NA ),
    axis.line = element_blank()
  )+
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = expression(
      italic("δ") * "(" * italic(Y) * ")" ==
        2.5 + 0.3 * italic(Y)
    ),
    parse = TRUE,
    hjust = -0.7,
    vjust = 1.7,
    size = 5,
    family = "Arial"
  )+
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = sprintf("MAE[rejected] == %.1f", SOM_MAE$MAE_rejected),
    parse = TRUE,
    hjust = 1.05,
    vjust = -2.2,
    colour = "#67000D",
    size = 5,
    family = "Arial"
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = sprintf("MAE[accepted] == %.1f", SOM_MAE$MAE_accepted),
    parse = TRUE,
    hjust = 1.05,
    vjust = -0.8,
    colour = "#1F4E79",
    size = 5,
    family = "Arial"
  ) +
  facet_grid(. ~ title, labeller = label_parsed)

SOM_Rejector_1to1







#----------------------------------------------------------------------------------------------------------------
# K











K_Meh3_threshold_alpha <- 0.05
K_Meh3_abs_buffer  <- 25
K_Meh3_rel_factor  <- 0.3




K_Meh3_p_value <- compute_p_exceed_hybrid(
  quantile_matrix = K_Meh3_pred_quantiles,
  pred_values     = K_Meh3_pred,
  a               = K_Meh3_abs_buffer,
  b               = K_Meh3_rel_factor,
  quantile_probs  = quantiles
)

K_Meh3_p_value


K_Meh3_error_p.df <- data.frame(
  pred    = K_Meh3_pred$K_Meh3_pred,
  true    = K_Meh3_true$K_Meh3_target,
  p_value = K_Meh3_p_value$p_value
) %>%
  mutate(
    pred  = as.numeric(pred),
    true  = as.numeric(true),
    error = abs(pred - true),

    # ---- thresholds ----
    delta_model = K_Meh3_abs_buffer + K_Meh3_rel_factor * pred,  # for model-side probability
    delta_true  = K_Meh3_abs_buffer + K_Meh3_rel_factor * true,  # REAL boundary for outliers

    # ---- evaluation ----
    is_wrong    = error > delta_true,
    is_rejected = p_value > K_Meh3_threshold_alpha
  )

K_Meh3_error_p.df



K_Meh3_confusion <- factor(
  ifelse(K_Meh3_error_p.df$is_rejected &  K_Meh3_error_p.df$is_wrong, "Correct Rejection",
         ifelse(K_Meh3_error_p.df$is_rejected & !K_Meh3_error_p.df$is_wrong, "Incorrect Rejection",
                ifelse(!K_Meh3_error_p.df$is_rejected &  K_Meh3_error_p.df$is_wrong, "Incorrect Acceptance",
                       "Correct Acceptance"))),
  levels = c(
    "Incorrect Rejection",
    "Incorrect Acceptance",
    "Correct Rejection",
    "Correct Acceptance"
  )
)



K_Meh3_error_cor_plot.df <- K_Meh3_error_p.df %>%
  mutate(
    confusion = K_Meh3_confusion,
    title = "bold(K~(mg~kg^{-1})*' - TabPFN')"
  )

confusion_labels_K_Meh3 <- K_Meh3_error_cor_plot.df %>%
  mutate(
    group = ifelse(is_rejected, "Rejected", "Accepted")
  ) %>%
  count(confusion, group) %>%
  group_by(group) %>%
  mutate(
    pct_within = n / sum(n) * 100
  ) %>%
  ungroup() %>%
  mutate(
    pct_total = n / sum(n) * 100
  ) %>%
  mutate(
    label = sprintf(
      "%s: %.1f%% (%.1f%% of %s)",
      confusion,
      pct_total,
      pct_within,
      tolower(group)
    )
  )

K_Meh3_Rejector_AE_alpha <- ggplot(
  K_Meh3_error_cor_plot.df,
  aes(x = p_value,
      y = error,
      colour = confusion)
) +
  geom_point(
    shape = 21,
    aes(fill = confusion),
    colour = "black",
    stroke = 0.4,
    alpha = 0.5,
    size = 2.5
  ) +
  geom_vline(
    xintercept = K_Meh3_threshold_alpha,
    linetype = "dashed"
  ) +

  scale_fill_manual(
    values = c(
      "Incorrect Rejection" = "#67000D",
      "Incorrect Acceptance" = "#F4A6A6",
      "Correct Rejection"  = "#6BAED6",
      "Correct Acceptance"  = "#1F4E79"
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_log10(
    breaks = c(0.03,0.3, 3,30,300),
    labels = c("0.03","0.3", "3", "30", "300"),
    expand = expansion(mult = c(0.02, 0.02))
  )+
  labs(
    y = "Absolute Error",
    x = expression(
      plain("Pr") *
        "(" *
        "|" * italic(Y) * " − " * italic("ŷ") * "|" *
        italic(" ≥ ") *
        italic("δ") * "(" * italic(Y) * ")" *
        "|" * italic(x) *
        ")")) +
  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 20, hjust = 1.3),
    axis.title.y = element_blank(),
    strip.background = element_rect(fill = "grey95",linewidth = 0.5 ),
    legend.position = "none",
    strip.text = element_text(size = 18, face = "bold")
  )+
  annotate(
    "text",
    x = K_Meh3_threshold_alpha + 0.003,
    y = Inf,
    label = bquote(alpha == .(K_Meh3_threshold_alpha)),
    vjust = 1.5,
    hjust = 0,
    size = 5,
    family = "Arial"
  ) +
  annotate(
    "text",
    x = Inf,
    y =  min(K_Meh3_error_cor_plot.df$error),
    label = confusion_labels_K_Meh3$label,
    hjust = 1.02,
    vjust = seq(-0.5, by = -1.2, length.out = nrow(confusion_labels_K_Meh3)),
    colour = c(
      "Incorrect Rejection" = "#67000D",
      "Incorrect Acceptance" = "#F4A6A6",
      "Correct Rejection"  = "#6BAED6",
      "Correct Acceptance"  = "#1F4E79"
    )[confusion_labels_K_Meh3$confusion],
    size = 5,
    family = "Arial"
  )+
  facet_grid(. ~ title, labeller = label_parsed)

K_Meh3_Rejector_AE_alpha




K_Meh3_MAE <- K_Meh3_error_cor_plot.df %>%
  summarise(
    MAE_accepted = mean(error[!is_rejected], na.rm = TRUE),
    MAE_rejected = mean(error[is_rejected], na.rm = TRUE)
  )


K_Meh3_error_cor_plot.df$title <- "bold(K~(mg~kg^{-1})*' - TabPFN')"


K_Meh3_Rejector_1to1 <- ggplot(
  K_Meh3_error_cor_plot.df,
  aes(x = true, y = pred, fill = p_value)
) +
  geom_point(
    shape = 21,
    colour = "black",
    stroke = 0.4,
    size = 2.5,
    alpha = 0.9
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 1
  ) +

  geom_abline(
    intercept = K_Meh3_abs_buffer,
    slope = 1 + K_Meh3_rel_factor,
    linetype = "dotted",
    linewidth = 0.9
  ) +
  geom_abline(
    intercept = -K_Meh3_abs_buffer,
    slope = 1 - K_Meh3_rel_factor,
    linetype = "dotted",
    linewidth = 0.9
  ) +

  scale_fill_gradientn(
    colours = c(
      "#1F4E79",
      "#6BAED6",
      "#B2E2E2",
      "#FEC44F",
      "#FB6A4A",
      "darkred",
      "#67000D"
    ),
    values = scales::rescale(
      c(0, 0.015, 0.05, 0.0501, 0.2, 0.4, 1),
      from = c(0, 1)
    ),
    limits = c(0, 1),
    breaks = c(0, 0.05, 1),
    labels = c(0, expression(alpha == 0.05), 1)
  ) +
  scale_x_continuous(limits = c(0, 500), expand = c(0.005, 0.005)) +
  scale_y_continuous(limits = c(0, 500), expand = c(0.005, 0.005)) +
  labs(
    x = "Measured Value",
    y = "Predicted Value"
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 20, hjust = 1.3),
    axis.title.y = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill = "grey95",linewidth = 0.5 ),
    strip.text = element_text(size = 18, face = "bold")
  )+
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = bquote(
      italic("δ")(italic(Y)) ==
        .(K_Meh3_abs_buffer) +
        .(format(K_Meh3_rel_factor, trim = TRUE)) * italic(Y)
    ),
    hjust = -0.7,
    vjust = 1.7,
    size = 5,
    family = "Arial"
  )+
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = sprintf("MAE[rejected] == %.1f", K_Meh3_MAE$MAE_rejected),
    parse = TRUE,
    hjust = 1.05,
    vjust = -2.2,
    colour = "#67000D",
    size = 5,
    family = "Arial"
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = sprintf("MAE[accepted] == %.1f", K_Meh3_MAE$MAE_accepted),
    parse = TRUE,
    hjust = 1.05,
    vjust = -0.8,
    colour = "#1F4E79",
    size = 5,
    family = "Arial"
  ) +
  facet_grid(. ~ title, labeller = label_parsed)

K_Meh3_Rejector_1to1












#----------------------------------------------------------------------------------------------------------------
# P











P_Meh3_threshold_alpha <- 0.05
P_Meh3_abs_buffer  <- 10
P_Meh3_rel_factor  <- 0.3




P_Meh3_p_value <- compute_p_exceed_hybrid(
  quantile_matrix = P_Meh3_pred_quantiles,
  pred_values     = P_Meh3_pred,
  a               = P_Meh3_abs_buffer,
  b               = P_Meh3_rel_factor,
  quantile_probs  = quantiles
)

P_Meh3_p_value


P_Meh3_error_p.df <- data.frame(
  pred    = P_Meh3_pred$P_Meh3_pred,
  true    = P_Meh3_true$P_Meh3_target,
  p_value = P_Meh3_p_value$p_value
) %>%
  mutate(
    pred  = as.numeric(pred),
    true  = as.numeric(true),
    error = abs(pred - true),

    delta_model = P_Meh3_abs_buffer + P_Meh3_rel_factor * pred,
    delta_true  = P_Meh3_abs_buffer + P_Meh3_rel_factor * true,

    is_wrong    = error > delta_true,
    is_rejected = p_value > P_Meh3_threshold_alpha
  )

P_Meh3_error_p.df



P_Meh3_confusion <- factor(
  ifelse(P_Meh3_error_p.df$is_rejected &  P_Meh3_error_p.df$is_wrong, "Correct Rejection",
         ifelse(P_Meh3_error_p.df$is_rejected & !P_Meh3_error_p.df$is_wrong, "Incorrect Rejection",
                ifelse(!P_Meh3_error_p.df$is_rejected &  P_Meh3_error_p.df$is_wrong, "Incorrect Acceptance",
                       "Correct Acceptance"))),
  levels = c(
    "Incorrect Rejection",
    "Incorrect Acceptance",
    "Correct Rejection",
    "Correct Acceptance"
  )
)



P_Meh3_error_cor_plot.df <- P_Meh3_error_p.df %>%
  mutate(
    confusion = P_Meh3_confusion,
    title = "bold(P~(mg~kg^{-1})*' - TabPFN')"
  )

confusion_labels_P_Meh3 <- P_Meh3_error_cor_plot.df %>%
  mutate(
    group = ifelse(is_rejected, "Rejected", "Accepted")
  ) %>%
  count(confusion, group) %>%
  group_by(group) %>%
  mutate(
    pct_within = n / sum(n) * 100
  ) %>%
  ungroup() %>%
  mutate(
    pct_total = n / sum(n) * 100
  ) %>%
  mutate(
    label = sprintf(
      "%s: %.1f%% (%.1f%% of %s)",
      confusion,
      pct_total,
      pct_within,
      tolower(group)
    )
  )

P_Meh3_error_cor_plot.df$title <- "bold(P~(mg~kg^{-1})*' - TabPFN')"

P_Meh3_Rejector_AE_alpha <- ggplot(
  P_Meh3_error_cor_plot.df,
  aes(x = p_value,
      y = error,
      colour = confusion)
) +
  geom_point(
    shape = 21,
    aes(fill = confusion),
    colour = "black",
    stroke = 0.4,
    alpha = 0.5,
    size = 2.5
  ) +
  geom_vline(
    xintercept = P_Meh3_threshold_alpha,
    linetype = "dashed"
  ) +

  scale_fill_manual(
    values = c(
      "Incorrect Rejection" = "#67000D",
      "Incorrect Acceptance" = "#F4A6A6",
      "Correct Rejection"  = "#6BAED6",
      "Correct Acceptance"  = "#1F4E79"
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_log10(
    breaks = c(0.03,0.3, 3,30,300),
    labels = c("0.03","0.3", "3", "30", "300"),
    expand = expansion(mult = c(0.02, 0.02))
  )+


  labs(
    y = "Absolute Error",
    x = expression(italic(P)(abs(italic(Y) - hat(italic(y))) >= delta))
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 14,color="black"),
    axis.title =  element_blank(),
    legend.text = element_text(size = 18,color="black"),
    legend.title = element_blank(),
    strip.text = element_text(size = 18,colour = "black", face = "bold"),
    strip.background = element_rect(fill = "grey95",linewidth = 0.5 ),
    panel.border = element_rect(colour = "black", linewidth = 0.5,fill = NA ),
    axis.line = element_blank(),
    legend.position = "none"
  )+
  annotate(
    "text",
    x = P_Meh3_threshold_alpha + 0.003,
    y = Inf,
    label = bquote(alpha == .(P_Meh3_threshold_alpha)),
    vjust = 1.5,
    hjust = 0,
    size = 5,
    family = "Arial"
  ) +
  annotate(
    "text",
    x = Inf,
    y =  min(P_Meh3_error_cor_plot.df$error),
    label = confusion_labels_P_Meh3$label,
    hjust = 1.02,
    vjust = seq(-0.5, by = -1.2, length.out = nrow(confusion_labels_P_Meh3)),
    colour = c(
      "Incorrect Rejection" = "#67000D",
      "Incorrect Acceptance" = "#F4A6A6",
      "Correct Rejection"  = "#6BAED6",
      "Correct Acceptance"  = "#1F4E79"
    )[confusion_labels_P_Meh3$confusion],
    size = 5,
    family = "Arial"
  )+
  facet_grid(. ~ title, labeller = label_parsed)

P_Meh3_Rejector_AE_alpha



P_Meh3_MAE <- P_Meh3_error_cor_plot.df %>%
  summarise(
    MAE_accepted = mean(error[!is_rejected], na.rm = TRUE),
    MAE_rejected = mean(error[is_rejected], na.rm = TRUE)
  )


P_Meh3_error_cor_plot.df$title <- "bold(P~(mg~kg^{-1})*' - TabPFN')"


P_Meh3_Rejector_1to1 <- ggplot(
  P_Meh3_error_cor_plot.df,
  aes(x = true, y = pred, fill = p_value)
) +
  geom_point(
    shape = 21,
    colour = "black",
    stroke = 0.4,
    size = 2.5,
    alpha = 0.9
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 1
  ) +

  geom_abline(
    intercept = P_Meh3_abs_buffer,
    slope = 1 + P_Meh3_rel_factor,
    linetype = "dotted",
    linewidth = 0.9
  ) +
  geom_abline(
    intercept = -P_Meh3_abs_buffer,
    slope = 1 - P_Meh3_rel_factor,
    linetype = "dotted",
    linewidth = 0.9
  ) +

  scale_fill_gradientn(
    colours = c(
      "#1F4E79",
      "#6BAED6",
      "#B2E2E2",
      "#FEC44F",
      "#FB6A4A",
      "darkred",
      "#67000D"
    ),
    values = scales::rescale(
      c(0, 0.015, 0.05, 0.0501, 0.2, 0.4, 1),
      from = c(0, 1)
    ),
    limits = c(0, 1),
    breaks = c(0, 0.05, 1),
    labels = c(0, expression(alpha == 0.05), 1)
  ) +
  scale_x_continuous(limits = c(0, 430),expand = c(0.005, 0.005)) +
  scale_y_continuous(limits = c(0, 430),expand = c(0.005, 0.005)) +
  labs(
    x = "Measured Value",
    y = "Predicted Value"
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 14,color="black"),
    axis.title =  element_blank(),
    legend.text = element_text(size = 14,color="black"),
    legend.title = element_text(size = 18,color="black"),
    strip.text = element_text(size = 18,colour = "black", face = "bold"),
    strip.background = element_rect(fill = "grey95",linewidth = 0.5 ),
    panel.border = element_rect(colour = "black", linewidth = 0.5,fill = NA ),
    axis.line = element_blank(),
    legend.position = "none"
  )+
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = bquote(
      italic("δ")(italic(Y)) ==
        .(P_Meh3_abs_buffer) +
        .(format(P_Meh3_rel_factor, trim = TRUE)) * italic(Y)
    ),
    hjust = -0.7,
    vjust = 1.7,
    size = 5,
    family = "Arial"
  )+
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = sprintf("MAE[rejected] == %.1f", P_Meh3_MAE$MAE_rejected),
    parse = TRUE,
    hjust = 1.05,
    vjust = -2.2,
    colour = "#67000D",
    size = 5,
    family = "Arial"
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = sprintf("MAE[accepted] == %.1f", P_Meh3_MAE$MAE_accepted),
    parse = TRUE,
    hjust = 1.05,
    vjust = -0.8,
    colour = "#1F4E79",
    size = 5,
    family = "Arial"
  ) +
  facet_grid(. ~ title, labeller = label_parsed)

P_Meh3_Rejector_1to1










Rejector_Matrix_AE_alpha_p_relative <- egg::ggarrange(
  Clay_Rejector_AE_alpha,
  SOM_Rejector_AE_alpha,
  K_Meh3_Rejector_AE_alpha,
  P_Meh3_Rejector_AE_alpha, ncol=2,nrow=2)


Rejector_Matrix_1to1_p_relative <- egg::ggarrange(
  Clay_Rejector_1to1,
  SOM_Rejector_1to1,
  K_Meh3_Rejector_1to1,
  P_Meh3_Rejector_1to1, ncol=2,nrow=2)

ggsave("Figures/Fig3.jpg", plot = Rejector_Matrix_AE_alpha_p_relative, width = 3000*1.5, height = 2600*1.5, units = "px", dpi = 300)
ggsave("Figures/Fig2.jpg", plot = Rejector_Matrix_1to1_p_relative, width = 3000*1.5, height = 2600*1.5, units = "px", dpi = 300)

