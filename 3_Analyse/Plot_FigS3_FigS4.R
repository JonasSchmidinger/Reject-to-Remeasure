# In this script, we create Fig.S3 and Fig.S4

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


Total_C_true<- read.csv(file = "Results/TabPFN/Outer/Total_C_true.csv")
Sand_true<- read.csv(file = "Results/TabPFN/Outer/Sand_true.csv")
pH_SMP_true<-  read.csv(file = "Results/TabPFN/Outer/pH_SMP_true.csv")
Al_Meh3_true<- read.csv(file = "Results/TabPFN/Outer/Al_Meh3_true.csv")


Total_C_pred<- read.csv(file = "Results/TabPFN/Outer/Total_C_pred.csv")
Sand_pred<- read.csv(file = "Results/TabICL/Outer/Sand_pred.csv")
pH_SMP_pred<-  read.csv(file = "Results/TabPFN/Outer/pH_SMP_pred.csv")
Al_Meh3_pred<- read.csv(file = "Results/TabICL/Outer/Al_Meh3_pred.csv")


Total_C_pred_quantiles<- read.csv(file = "Results/TabPFN/Outer/Total_C_pred_quantiles.csv")
Sand_pred_quantiles<- read.csv(file = "Results/TabICL/Outer/Sand_pred_quantiles.csv")
pH_SMP_pred_quantiles<-  read.csv(file = "Results/TabPFN/Outer/pH_SMP_pred_quantiles.csv")
Al_Meh3_pred_quantiles<- read.csv(file = "Results/TabICL/Outer/Al_Meh3_pred_quantiles.csv")




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

    #Need different criteria not based on y_hat but on Y, as changes along Y will also change the threshold (as compared to absolute threshold with fixed threshold widths)
    L <- (y_hat - a) / (1 + b)
    U <- (y_hat + a) / (1 - b)

    p_vals[i] <- cdf_fun(L) + (1 - cdf_fun(U))
  }

  data.frame(p_value = p_vals)
}



#----------------------------------------------------------------------------------------------------------------
# Sand


Sand_abs_buffer  <- 120
Sand_rel_factor  <- 0
Sand_threshold_alpha <- 0.05






Sand_p_value <- compute_p_exceed_hybrid(
  quantile_matrix = Sand_pred_quantiles * 10,
  pred_values     = data.frame(pred = Sand_pred$Sand_pred * 10),
  a               = Sand_abs_buffer,
  b               = Sand_rel_factor,
  quantile_probs  = quantiles
)







Sand_error_p.df <- data.frame(
  pred    =  as.numeric(Sand_pred$Sand_pred) * 10,
  true    = as.numeric(Sand_true$Sand_target) * 10,
  p_value = Sand_p_value$p_value
) %>%
  mutate(
    pred  = as.numeric(pred),
    true  = as.numeric(true),
    error = abs(pred - true),

    delta_model = Sand_abs_buffer + Sand_rel_factor * pred,
    delta_true  = Sand_abs_buffer + Sand_rel_factor * true,

    is_wrong    = error > delta_true,
    is_rejected = p_value > Sand_threshold_alpha
  )

Sand_error_p.df



Sand_confusion <- factor(
  ifelse(Sand_error_p.df$is_rejected &  Sand_error_p.df$is_wrong, "Correct Rejection",
         ifelse(Sand_error_p.df$is_rejected & !Sand_error_p.df$is_wrong, "Incorrect Rejection",
                ifelse(!Sand_error_p.df$is_rejected &  Sand_error_p.df$is_wrong, "Incorrect Acceptance",
                       "Correct Acceptance"))),
  levels = c(
    "Incorrect Rejection",
    "Incorrect Acceptance",
    "Correct Rejection",
    "Correct Acceptance"
  )
)

Sand_error_cor_plot.df <- Sand_error_p.df %>%
  mutate(
    confusion = Sand_confusion,
    title =  "bold(Sand~(g~kg^{-1})*' - TabICL')"
  )

confusion_labels_Sand <- Sand_error_cor_plot.df %>%
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
confusion_labels_Sand

acceptance_rate <- Sand_error_cor_plot.df %>%
  summarise(rate = mean(!is_rejected) * 100) %>%
  pull(rate)

rejection_rate <- Sand_error_cor_plot.df %>%
  summarise(rate = mean(is_rejected) * 100) %>%
  pull(rate)


Sand_error_cor_plot.df$p_value

Sand_Rejector_AE_alpha <-ggplot(
  Sand_error_cor_plot.df,
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
    label = confusion_labels_Sand$label,
    hjust = 1.02,
    vjust = seq(1.2, by = 1.3, length.out = nrow(confusion_labels_Sand)),
    colour = c(
      "Incorrect Rejection" = "#67000D",
      "Incorrect Acceptance" = "#F4A6A6",
      "Correct Rejection"  = "#6BAED6",
      "Correct Acceptance"  = "#1F4E79"
    )[confusion_labels_Sand$confusion],
    size = 5,
    family = "Arial"
  )+
  facet_grid(. ~ title, labeller = label_parsed)

Sand_Rejector_AE_alpha



Sand_MAE <- Sand_error_cor_plot.df %>%
  summarise(
    MAE_accepted = mean(error[!is_rejected], na.rm = TRUE),
    MAE_rejected = mean(error[is_rejected], na.rm = TRUE)
  )



Sand_error_cor_plot.df$title <- "bold(Sand~(g~kg^{-1})*' - TabICL')"

Sand_Rejector_1to1 <- ggplot(
  Sand_error_cor_plot.df,
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
  fun = function(x) x + (Sand_abs_buffer + Sand_rel_factor * x),
  linetype = "dotted",
  linewidth = 0.9
) +
  geom_function(
    fun = function(x) x - (Sand_abs_buffer + Sand_rel_factor * x),
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
        120 * italic(Y)
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
    label = sprintf("MAE[rejected] == %.1f", Sand_MAE$MAE_rejected),
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
    label = sprintf("MAE[accepted] == %.1f", Sand_MAE$MAE_accepted),
    parse = TRUE,
    hjust = 1.05,
    vjust = -0.8,
    colour = "#1F4E79",
    size = 5,
    family = "Arial"
  ) +
  facet_grid(. ~ title, labeller = label_parsed)

Sand_Rejector_1to1





#----------------------------------------------------------------------------------------------------------------
# Total_C




Total_C_abs_buffer  <- 10
Total_C_rel_factor  <- 0
Total_C_threshold_alpha <- 0.05



Total_C_p_value <- compute_p_exceed_hybrid(
  quantile_matrix = Total_C_pred_quantiles * 10,
  pred_values     = data.frame(pred = Total_C_pred$Total_C_pred * 10),
  a               = Total_C_abs_buffer,
  b               = Total_C_rel_factor,
  quantile_probs  = quantiles
)







Total_C_error_p.df <- data.frame(
  pred    =  as.numeric(Total_C_pred$Total_C_pred) * 10,
  true    = as.numeric(Total_C_true$Total_C_target) * 10,
  p_value = Total_C_p_value$p_value
) %>%
  mutate(
    pred  = as.numeric(pred),
    true  = as.numeric(true),
    error = abs(pred - true),

    delta_model = Total_C_abs_buffer + Total_C_rel_factor * pred,
    delta_true  = Total_C_abs_buffer + Total_C_rel_factor * true,

    is_wrong    = error > delta_true,
    is_rejected = p_value > Total_C_threshold_alpha
  )

Total_C_error_p.df



Total_C_confusion <- factor(
  ifelse(Total_C_error_p.df$is_rejected &  Total_C_error_p.df$is_wrong, "Correct Rejection",
         ifelse(Total_C_error_p.df$is_rejected & !Total_C_error_p.df$is_wrong, "Incorrect Rejection",
                ifelse(!Total_C_error_p.df$is_rejected &  Total_C_error_p.df$is_wrong, "Incorrect Acceptance",
                       "Correct Acceptance"))),
  levels = c(
    "Incorrect Rejection",
    "Incorrect Acceptance",
    "Correct Rejection",
    "Correct Acceptance"
  )
)



Total_C_error_cor_plot.df <- Total_C_error_p.df %>%
  mutate(
    confusion = Total_C_confusion,
    title = "bold(TC~(g~kg^{-1})*' - TabPFN')"
  )

confusion_labels_Total_C <- Total_C_error_cor_plot.df %>%
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
confusion_labels_Total_C




Total_C_error_cor_plot.df$p_value


Total_C_Rejector_AE_alpha <-ggplot(
  Total_C_error_cor_plot.df,
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
    label = confusion_labels_Total_C$label,
    hjust = 1.02,
    vjust = seq(1.2, by = 1.3, length.out = nrow(confusion_labels_Total_C)),
    colour = c(
      "Incorrect Rejection" = "#67000D",
      "Incorrect Acceptance" = "#F4A6A6",
      "Correct Rejection"  = "#6BAED6",
      "Correct Acceptance"  = "#1F4E79"
    )[confusion_labels_Total_C$confusion],
    size = 5,
    family = "Arial"
  )+
  facet_grid(. ~ title, labeller = label_parsed)

Total_C_Rejector_AE_alpha



Total_C_MAE <- Total_C_error_cor_plot.df %>%
  summarise(
    MAE_accepted = mean(error[!is_rejected], na.rm = TRUE),
    MAE_rejected = mean(error[is_rejected], na.rm = TRUE)
  )

Total_C_error_cor_plot.df$title <- "bold(TC~(g~kg^{-1})*' - TabPFN')"
Total_C_error_cor_plot.df

Total_C_Rejector_1to1 <- ggplot(
  Total_C_error_cor_plot.df,
  aes(x = true, y = pred, fill = p_value)) +
  geom_point(
    shape = 21,
    colour = "black",
    stroke = 0.4,
    size = 2.5,
    alpha = 0.9
  )+
  # 1:1 line
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 1
  ) +

  geom_function(
    fun = function(x) x + (Total_C_abs_buffer + Total_C_rel_factor * x),
    linetype = "dotted",
    linewidth = 0.9
  ) +
  geom_function(
    fun = function(x) x - (Total_C_abs_buffer + Total_C_rel_factor * x),
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
  scale_x_continuous(limits = c(0, 120),expand = c(0.005, 0.005)) +
  scale_y_continuous(limits = c(0, 120),expand = c(0.005, 0.005)) +

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
        10
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
    label = sprintf("MAE[rejected] == %.1f", Total_C_MAE$MAE_rejected),
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
    label = sprintf("MAE[accepted] == %.1f", Total_C_MAE$MAE_accepted),
    parse = TRUE,
    hjust = 1.05,
    vjust = -0.8,
    colour = "#1F4E79",
    size = 5,
    family = "Arial"
  ) +
  facet_grid(. ~ title, labeller = label_parsed)

Total_C_Rejector_1to1



#----------------------------------------------------------------------------------------------------------------
# Al

Al_Meh3_threshold_alpha <- 0.05
Al_Meh3_abs_buffer  <- 250
Al_Meh3_rel_factor  <- 0




Al_Meh3_p_value <- compute_p_exceed_hybrid(
  quantile_matrix = Al_Meh3_pred_quantiles,
  pred_values     = Al_Meh3_pred,
  a               = Al_Meh3_abs_buffer,
  b               = Al_Meh3_rel_factor,
  quantile_probs  = quantiles
)

Al_Meh3_p_value


Al_Meh3_error_p.df <- data.frame(
  pred    = Al_Meh3_pred$Al_Meh3_pred,
  true    = Al_Meh3_true$Al_Meh3_target,
  p_value = Al_Meh3_p_value$p_value
) %>%
  mutate(
    pred  = as.numeric(pred),
    true  = as.numeric(true),
    error = abs(pred - true),

    delta_model = Al_Meh3_abs_buffer + Al_Meh3_rel_factor * pred,
    delta_true  = Al_Meh3_abs_buffer + Al_Meh3_rel_factor * true,

    is_wrong    = error > delta_true,
    is_rejected = p_value > Al_Meh3_threshold_alpha
  )

Al_Meh3_error_p.df



Al_Meh3_confusion <- factor(
  ifelse(Al_Meh3_error_p.df$is_rejected &  Al_Meh3_error_p.df$is_wrong, "Correct Rejection",
         ifelse(Al_Meh3_error_p.df$is_rejected & !Al_Meh3_error_p.df$is_wrong, "Incorrect Rejection",
                ifelse(!Al_Meh3_error_p.df$is_rejected &  Al_Meh3_error_p.df$is_wrong, "Incorrect Acceptance",
                       "Correct Acceptance"))),
  levels = c(
    "Incorrect Rejection",
    "Incorrect Acceptance",
    "Correct Rejection",
    "Correct Acceptance"
  )
)



Al_Meh3_error_cor_plot.df <- Al_Meh3_error_p.df %>%
  mutate(
    confusion = Al_Meh3_confusion,
    title = "bold(Al~(mg~kg^{-1})*' - TabICL')"
  )

confusion_labels_Al_Meh3 <- Al_Meh3_error_cor_plot.df %>%
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
Al_Meh3_error_cor_plot.df$title <- "bold(Al~(mg~kg^{-1})*' - TabICL')"

Al_Meh3_Rejector_AE_alpha <- ggplot(
  Al_Meh3_error_cor_plot.df,
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
    xintercept = Al_Meh3_threshold_alpha,
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
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
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
    x = Al_Meh3_threshold_alpha + 0.003,
    y = Inf,
    label = bquote(alpha == .(Al_Meh3_threshold_alpha)),
    vjust = 1.5,
    hjust = 0,
    size = 5,
    family = "Arial"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = confusion_labels_Al_Meh3$label,
    hjust = 1.02,
    vjust = seq(1.2, by = 1.3, length.out = nrow(confusion_labels_Al_Meh3)),
    colour = c(
      "Incorrect Rejection" = "#67000D",
      "Incorrect Acceptance" = "#F4A6A6",
      "Correct Rejection"  = "#6BAED6",
      "Correct Acceptance"  = "#1F4E79"
    )[confusion_labels_Al_Meh3$confusion],
    size = 5,
    family = "Arial"
  ) +
  facet_grid(. ~ title, labeller = label_parsed)

Al_Meh3_Rejector_AE_alpha


Al_Meh3_MAE <- Al_Meh3_error_cor_plot.df %>%
  summarise(
    MAE_accepted = mean(error[!is_rejected], na.rm = TRUE),
    MAE_rejected = mean(error[is_rejected], na.rm = TRUE)
  )



Al_Meh3_error_cor_plot.df$title <- "bold(Al~(mg~kg^{-1})*' - TabICL')"


Al_Meh3_Rejector_1to1 <- ggplot(
  Al_Meh3_error_cor_plot.df,
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
  intercept = Al_Meh3_abs_buffer,
  slope = 1 + Al_Meh3_rel_factor,
  linetype = "dotted",
  linewidth = 0.9
) +
  geom_abline(
    intercept = -Al_Meh3_abs_buffer,
    slope = 1 - Al_Meh3_rel_factor,
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
  scale_x_continuous(limits = c(0, 2500), expand = c(0.005, 0.005)) +
  scale_y_continuous(limits = c(0, 2500), expand = c(0.005, 0.005)) +
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
        .(Al_Meh3_abs_buffer)
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
    label = sprintf("MAE[rejected] == %.1f", Al_Meh3_MAE$MAE_rejected),
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
    label = sprintf("MAE[accepted] == %.1f", Al_Meh3_MAE$MAE_accepted),
    parse = TRUE,
    hjust = 1.05,
    vjust = -0.8,
    colour = "#1F4E79",
    size = 5,
    family = "Arial"
  ) +
  facet_grid(. ~ title, labeller = label_parsed)


Al_Meh3_Rejector_1to1


#----------------------------------------------------------------------------------------------------------------
# pH




pH_SMP_threshold = 0.25
pH_SMP_threshold_alpha = 0.05

pH_SMP_threshold_alpha <- 0.05
pH_SMP_abs_buffer      <-  pH_SMP_threshold



pH_SMP_p_value <- compute_p_exceed_all(
  quantile_matrix =  pH_SMP_pred_quantiles,
  pred_values     =  pH_SMP_pred,
  threshold       =  pH_SMP_abs_buffer,
  quantile_probs  = quantiles
)


pH_SMP_error_p.df <- data.frame(
  pred    =  pH_SMP_pred$ pH_SMP_pred,
  true    =  pH_SMP_true$ pH_SMP_target,
  p_value =  pH_SMP_p_value$p_value
) %>%
  mutate(
    pred  = as.numeric(pred),
    true  = as.numeric(true),
    error = abs(pred - true),

    delta_model =  pH_SMP_abs_buffer,
    delta_true  =  pH_SMP_abs_buffer,

    is_wrong    = error > delta_true,
    is_rejected = p_value >  pH_SMP_threshold_alpha
  )

pH_SMP_error_p.df


pH_SMP_confusion <- factor(
  ifelse( pH_SMP_error_p.df$is_rejected &   pH_SMP_error_p.df$is_wrong, "Correct Rejection",
          ifelse( pH_SMP_error_p.df$is_rejected & ! pH_SMP_error_p.df$is_wrong, "Incorrect Rejection",
                  ifelse(! pH_SMP_error_p.df$is_rejected &   pH_SMP_error_p.df$is_wrong, "Incorrect Acceptance",
                         "Correct Acceptance"))),
  levels = c(
    "Incorrect Rejection",
    "Incorrect Acceptance",
    "Correct Rejection",
    "Correct Acceptance"
  )
)


pH_SMP_error_cor_plot.df <-  pH_SMP_error_p.df %>%
  mutate(
    confusion =  pH_SMP_confusion,
    title = "bold(pH~' - TabPFN')"
  )

pH_SMP_error_cor_plot.df


confusion_labels_pH_SMP <- pH_SMP_error_cor_plot.df %>%
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
confusion_labels_pH_SMP


pH_SMP_Rejector_AE_alpha <- ggplot(
  pH_SMP_error_cor_plot.df,
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
    xintercept = pH_SMP_threshold_alpha,
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
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
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
    x = pH_SMP_threshold_alpha + 0.003,
    y = Inf,
    label = bquote(alpha == .(pH_SMP_threshold_alpha)),
    vjust = 1.5,
    hjust = 0,
    size = 5,
    family = "Arial"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = confusion_labels_pH_SMP$label,
    hjust = 1.02,
    vjust = seq(1.2, by = 1.3, length.out = nrow(confusion_labels_pH_SMP)),
    colour = c(
      "Incorrect Rejection" = "#67000D",
      "Incorrect Acceptance" = "#F4A6A6",
      "Correct Rejection"  = "#6BAED6",
      "Correct Acceptance"  = "#1F4E79"
    )[confusion_labels_pH_SMP$confusion],
    size = 5,
    family = "Arial"
  ) +
  facet_grid(. ~ title, labeller = label_parsed)

pH_SMP_Rejector_AE_alpha



pH_SMP_MAE <- pH_SMP_error_cor_plot.df %>%
  summarise(
    MAE_accepted = mean(error[!is_rejected], na.rm = TRUE),
    MAE_rejected = mean(error[is_rejected], na.rm = TRUE)
  )



pH_SMP_error_cor_plot.df$title <- "bold(pH~' - TabPFN')"


pH_SMP_Rejector_1to1 <- ggplot(
  pH_SMP_error_cor_plot.df,
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
  intercept =  pH_SMP_threshold,
  slope = 1,
  linetype = "dotted",
  linewidth = 0.9
) +
  geom_abline(
    intercept = - pH_SMP_threshold,
    slope = 1,
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
  scale_x_continuous(limits = c(4.8, 7.8),expand = c(0.005, 0.005)) +
  scale_y_continuous(limits = c(4.8, 7.8),expand = c(0.005, 0.005)) +
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
      italic("δ") * "(" * italic(Y) * ")" == .( pH_SMP_threshold)
    ),
    hjust = -0.7,
    vjust = 1.7,
    size = 5,
    family = "Arial"
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = sprintf("MAE[rejected] == %.2f", pH_SMP_MAE$MAE_rejected),
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
    label = sprintf("MAE[accepted] == %.2f", pH_SMP_MAE$MAE_accepted),
    parse = TRUE,
    hjust = 1.05,
    vjust = -0.8,
    colour = "#1F4E79",
    size = 5,
    family = "Arial"
  ) +
  facet_grid(. ~ title, labeller = label_parsed)

pH_SMP_Rejector_1to1













Rejector_Matrix_AE_alpha_p_relative_appendix <- egg::ggarrange(
  Sand_Rejector_AE_alpha,
  Total_C_Rejector_AE_alpha,
  Al_Meh3_Rejector_AE_alpha,
  pH_SMP_Rejector_AE_alpha, ncol=2,nrow=2)

Rejector_Matrix_AE_alpha_p_relative_appendix

Rejector_Matrix_1to1_p_relative_appendix  <- egg::ggarrange(
  Sand_Rejector_1to1,
  Total_C_Rejector_1to1,
  Al_Meh3_Rejector_1to1,
  pH_SMP_Rejector_1to1, ncol=2,nrow=2)

Rejector_Matrix_1to1_p_relative_appendix

ggsave("Figures/FigS4.jpg", plot = Rejector_Matrix_AE_alpha_p_relative_appendix, width = 3000*1.5, height = 2600*1.5, units = "px", dpi = 300)
ggsave("Figures/FigS3.jpg", plot = Rejector_Matrix_1to1_p_relative_appendix, width = 3000*1.5, height = 2600*1.5, units = "px", dpi = 300)

