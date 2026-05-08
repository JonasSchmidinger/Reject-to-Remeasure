# In this script, we create Fig.4

library(tidyverse)
library(ggplot2)
library(dplyr)
library(scales)
library(future.apply)
library(showtext)
library(metR)

plan(multisession, workers = 6)


font_add(
  family = "Arial",
  regular = "C:/Windows/Fonts/arial.ttf",
  italic  = "C:/Windows/Fonts/ariali.ttf",
  bold    = "C:/Windows/Fonts/arialbd.ttf",
  bolditalic = "C:/Windows/Fonts/arialbi.ttf"
)

showtext_auto()


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




#----------------------------------------------------------------------------------------------------------------
# Clay

Clay_abs_buffer      <- 50
Clay_rel_factor      <- 0.15
Clay_threshold_alpha <- 0.05


Clay_p_value <- compute_p_exceed_hybrid(
  quantile_matrix = as.matrix(Clay_pred_quantiles) * 10,
  pred_values     = Clay_pred$Clay_pred * 10,
  a               = Clay_abs_buffer,
  b               = Clay_rel_factor,
  quantile_probs  = quantiles
)


Clay_error_p.df <- data.frame(
  pred    = as.numeric(Clay_pred$Clay_pred) * 10,
  true    = as.numeric(Clay_true$Clay_target) * 10,
  p_value = Clay_p_value$p_value
) %>%
  mutate(
    error = abs(pred - true),

    delta_true = Clay_abs_buffer + Clay_rel_factor * true,

    is_wrong = error > delta_true,

    is_rejected = p_value > Clay_threshold_alpha
  )



Clay_acceptance_rate <- mean(!Clay_error_p.df$is_rejected)



visnir_costs  <- seq(0.75, 10.7, by = 0.025)
texture_costs <- seq(33, 56.5, by = 0.025) # Slightly adjusted range for better visibility (done for all ranges !)

Clay_economic_cost <- expand.grid(
  visnir_cost  = visnir_costs,
  texture_cost = texture_costs
)


Clay_economic_cost <- Clay_economic_cost %>%
  mutate(

    baseline_cost = texture_cost,

    strategy_cost = visnir_cost +
      (1 - Clay_acceptance_rate) * texture_cost,

    cost_difference = strategy_cost - baseline_cost
  )


Clay_acceptance_percent <- round(Clay_acceptance_rate * 100, 1)

Clay_economic_cost$title <- paste0(
  "bold(Clay~'-'~'",
  Clay_acceptance_percent,
  "% Acceptance Rate')"
)



Clay_economic_cost <- Clay_economic_cost %>%
  mutate(
    visnir_cost  = as.numeric(visnir_cost),
    texture_cost = as.numeric(texture_cost),
    cost_difference  = as.numeric(cost_difference)
  )


contour_levels <- seq(-40, 40, by = 5)


Clay_economic_cost_plot <- ggplot(
  Clay_economic_cost,
  aes(
    x = visnir_cost,
    y = texture_cost,
    fill = cost_difference
  ))+
  geom_raster() +

  geom_contour(
    aes(z = cost_difference),
    breaks = contour_levels,
    colour = "black",
    linewidth = 0.4
  ) +

  geom_text_contour(
    aes(z = cost_difference),
    breaks = contour_levels,
    size = 6,
    stroke = 0.1,
    family = "Arial",
    skip = 0,
    rotate = TRUE
  )+

  scale_fill_gradientn(
    colours = c(
      "#01665E",
      "#35978F",
      "#80CDC1",
      "#F7F7F7",
      "#DFC27D",
      "#BF812D",
      "#8C510A"
    ),
    values = scales::rescale(
      c(-40,-20, -1, 0,1,20 ,40)
    ),
    limits = c(-40, 40),
    name = "Cost saving\n(CAD sample^{-1})"
  )+


  scale_y_continuous(
      breaks = seq(35, max(Clay_economic_cost$texture_cost, na.rm = TRUE), by = 10),
    expand = c(0,0)
  )+

  scale_x_continuous(
    expand = c(0,0)
  ) +

  labs(
    y = expression("Conventional Laboratory Cost " ~ "(" * "$CAD" * " Sample"^{-1} * ")"))+

  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.title.y = element_text(size = 20, hjust = 3,vjust= 3),
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

Clay_economic_cost_plot




#----------------------------------------------------------------------------------------------------------------
# SOM

SOM_abs_buffer      <- 3
SOM_rel_factor      <- 0.3
SOM_threshold_alpha <- 0.05



SOM_p_value <- compute_p_exceed_hybrid(
  quantile_matrix = as.matrix(SOM_pred_quantiles) * 10,
  pred_values     = SOM_pred$SOM_pred * 10,
  a               = SOM_abs_buffer,
  b               = SOM_rel_factor,
  quantile_probs  = quantiles
)



SOM_error_p.df <- data.frame(
  pred    = as.numeric(SOM_pred$SOM_pred) * 10,
  true    = as.numeric(SOM_true$SOM_target) * 10,
  p_value = SOM_p_value$p_value
) %>%
  mutate(
    error = abs(pred - true),

    delta_true = SOM_abs_buffer + SOM_rel_factor * true,

    is_wrong = error > delta_true,

    is_rejected = p_value > SOM_threshold_alpha
  )



SOM_acceptance_rate <- mean(!SOM_error_p.df$is_rejected)



visnir_costs  <- seq(0.75, 10.7, by = 0.025)
SOM_costs <- seq(13, 26, by = 0.025)

SOM_economic_cost <- expand.grid(
  visnir_cost  = visnir_costs,
  SOM_cost = SOM_costs
)


SOM_economic_cost <- SOM_economic_cost %>%
  mutate(

    baseline_cost = SOM_cost,

    strategy_cost = visnir_cost +
      (1 - SOM_acceptance_rate) * SOM_cost,

    cost_difference = strategy_cost - baseline_cost
  )


SOM_acceptance_percent <- round(SOM_acceptance_rate * 100, 1)

SOM_economic_cost$title <- paste0(
  "bold(SOM~'-'~'",
  SOM_acceptance_percent,
  "% Acceptance Rate')"
)



SOM_economic_cost <- SOM_economic_cost %>%
  mutate(
    visnir_cost  = as.numeric(visnir_cost),
    SOM_cost = as.numeric(SOM_cost),
    cost_difference  = as.numeric(cost_difference)
  )


SOM_economic_cost_plot <-ggplot(
  SOM_economic_cost,
  aes(
    x = visnir_cost,
    y = SOM_cost,
    fill = cost_difference
  ))+
  geom_raster() +

  geom_contour(
    aes(z = cost_difference),
    breaks = contour_levels,
    colour = "black",
    linewidth = 0.4
  ) +

  geom_text_contour(
    aes(z = cost_difference),
    breaks = contour_levels,
    size = 6,
    stroke = 0.1,
    family = "Arial",
    skip = 0,
    rotate = TRUE
  )+

  scale_fill_gradientn(
    colours = c(
      "#01665E",
      "#35978F",
      "#80CDC1",
      "#F7F7F7",
      "#DFC27D",
      "#BF812D",
      "#8C510A"
    ),
    values = scales::rescale(
      c(-40,-20, -1, 0,1,20 ,40)
    ),
    limits = c(-40, 40),
    name = expression("Cost Difference: Reject-to-Remeasure vs. Conventional" ~ "(" * "$" * "CAD Sample"^{-1} * ")"))+
  scale_y_continuous(
    breaks = seq(15, max(SOM_economic_cost$SOM_cost, na.rm = TRUE), by = 5),
    expand = c(0, 0)
  ) +

  scale_x_continuous(
    expand = c(0, 0)
  ) +

  labs(
    x = expression(delta[0]),
    y = expression(italic(b))
  ) +

  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_blank(),
    legend.position = c(-0.63, 1.07),
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




SOM_economic_cost_plot




#----------------------------------------------------------------------------------------------------------------
# K



K_Meh3_abs_buffer      <- 25
K_Meh3_rel_factor      <- 0.3
K_Meh3_threshold_alpha <- 0.05


K_Meh3_p_value <- compute_p_exceed_hybrid(
  quantile_matrix = as.matrix(K_Meh3_pred_quantiles),
  pred_values     = K_Meh3_pred$K_Meh3_pred,
  a               = K_Meh3_abs_buffer,
  b               = K_Meh3_rel_factor,
  quantile_probs  = quantiles
)


K_Meh3_error_p.df <- data.frame(
  pred    = as.numeric(K_Meh3_pred$K_Meh3_pred),
  true    = as.numeric(K_Meh3_true$K_Meh3_target),
  p_value = K_Meh3_p_value$p_value
) %>%
  mutate(
    error = abs(pred - true),

    delta_true = K_Meh3_abs_buffer + K_Meh3_rel_factor * true,

    is_wrong = error > delta_true,

    is_rejected = p_value > K_Meh3_threshold_alpha
  )



K_Meh3_acceptance_rate <- mean(!K_Meh3_error_p.df$is_rejected)



visnir_costs  <- seq(0.75, 10.7, by = 0.025)
K_Meh3_costs <-  seq(10, 22, by = 0.025)

K_Meh3_economic_cost <- expand.grid(
  visnir_cost  = visnir_costs,
  K_Meh3_cost = K_Meh3_costs
)


K_Meh3_economic_cost <- K_Meh3_economic_cost %>%
  mutate(

    baseline_cost = K_Meh3_cost,

    strategy_cost = visnir_cost +
      (1 - K_Meh3_acceptance_rate) * K_Meh3_cost,

    cost_difference = strategy_cost - baseline_cost
  )


K_Meh3_acceptance_percent <- round(K_Meh3_acceptance_rate * 100, 1)

K_Meh3_economic_cost$title <- paste0(
  "bold(K~'-'~'",
  K_Meh3_acceptance_percent,
  "% Acceptance Rate')"
)



K_Meh3_economic_cost <- K_Meh3_economic_cost %>%
  mutate(
    visnir_cost  = as.numeric(visnir_cost),
    K_Meh3_cost = as.numeric(K_Meh3_cost),
    cost_difference  = as.numeric(cost_difference)
  )


K_economic_cost_plot <-ggplot(
  K_Meh3_economic_cost,
  aes(
    x = visnir_cost,
    y = K_Meh3_cost,
    fill = cost_difference
  ))+
  geom_raster() +

  geom_contour(
    aes(z = cost_difference),
    breaks = contour_levels,
    colour = "black",
    linewidth = 0.4
  ) +

  geom_text_contour(
    aes(z = cost_difference),
    breaks = contour_levels,
    size = 6,
    stroke = 0.1,
    family = "Arial",
    skip = 0,
    rotate = TRUE
  )+
  scale_fill_gradientn(
    colours = c(
      "#01665E",
      "#35978F",
      "#80CDC1",
      "#F7F7F7",
      "#DFC27D",
      "#BF812D",
      "#8C510A"
    ),
    values = scales::rescale(
      c(-40,-20, -1, 0,1,20 ,40)
    ),
    limits = c(-40, 40),
    name = expression("Cost Saving with VNIRS" ~ "(" * "$" * "CAD Sample"^{-1} * ")"))+

  scale_y_continuous(
    breaks = seq(10, max(K_Meh3_economic_cost$K_Meh3_cost, na.rm = TRUE), by = 5),
    expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +

  labs(
    x = expression("VNIRS Cost " ~ "(" * "$CAD" * " Sample"^{-1} * ")"))+
  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.title.x = element_text(size = 20,hjust = 1.8, face = "italic"),
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



K_economic_cost_plot





#----------------------------------------------------------------------------------------------------------------
# P


P_Meh3_abs_buffer      <- 10
P_Meh3_rel_factor      <- 0.3
P_Meh3_threshold_alpha <- 0.05



P_Meh3_p_value <- compute_p_exceed_hybrid(
  quantile_matrix = as.matrix(P_Meh3_pred_quantiles),
  pred_values     = P_Meh3_pred$P_Meh3_pred,
  a               = P_Meh3_abs_buffer,
  b               = P_Meh3_rel_factor,
  quantile_probs  = quantiles
)



P_Meh3_error_p.df <- data.frame(
  pred    = as.numeric(P_Meh3_pred$P_Meh3_pred),
  true    = as.numeric(P_Meh3_true$P_Meh3_target),
  p_value = P_Meh3_p_value$p_value
) %>%
  mutate(
    error = abs(pred - true),

    delta_true = P_Meh3_abs_buffer + P_Meh3_rel_factor * true,

    is_wrong = error > delta_true,

    is_rejected = p_value > P_Meh3_threshold_alpha
  )



P_Meh3_acceptance_rate <- mean(!P_Meh3_error_p.df$is_rejected)



visnir_costs  <- seq(0.75, 10.7, by = 0.025)
P_Meh3_costs <- seq(10, 22, by = 0.025)

P_Meh3_economic_cost <- expand.grid(
  visnir_cost  = visnir_costs,
  P_Meh3_cost = P_Meh3_costs
)


P_Meh3_economic_cost <- P_Meh3_economic_cost %>%
  mutate(

    baseline_cost = P_Meh3_cost,

    strategy_cost = visnir_cost +
      (1 - P_Meh3_acceptance_rate) * P_Meh3_cost,

    cost_difference = strategy_cost - baseline_cost
  )


P_Meh3_acceptance_percent <- round(P_Meh3_acceptance_rate * 100, 1)

P_Meh3_economic_cost$title <- paste0(
  "bold(P~'-'~'",
  P_Meh3_acceptance_percent,
  "% Acceptance Rate')"
)



P_Meh3_economic_cost <- P_Meh3_economic_cost %>%
  mutate(
    visnir_cost  = as.numeric(visnir_cost),
    P_Meh3_cost = as.numeric(P_Meh3_cost),
    cost_difference  = as.numeric(cost_difference)
  )


P_economic_cost_plot <-ggplot(
  P_Meh3_economic_cost,
  aes(
    x = visnir_cost,
    y = P_Meh3_cost,
    fill = cost_difference
  ))+
  geom_raster() +

  geom_contour(
    aes(z = cost_difference),
    breaks = contour_levels,
    colour = "black",
    linewidth = 0.4
  ) +

  geom_text_contour(
    aes(z = cost_difference),
    breaks = contour_levels,
    size = 6,
    stroke = 0.1,
    family = "Arial",
    skip = 0,
    rotate = TRUE
  )+
  scale_fill_gradientn(
    colours = c(
      "#01665E",
      "#35978F",
      "#80CDC1",
      "#F7F7F7",
      "#DFC27D",
      "#BF812D",
      "#8C510A"
    ),
    values = scales::rescale(
      c(-40,-20, -1, 0,1,20 ,40)
    ),
    limits = c(-40, 40),
    name = expression("Cost saving with VNIRS" ~ "(" * "$" * "CAD Measurment"^{-1} * ")"))+

  scale_y_continuous(
    breaks = seq(10, max(P_Meh3_economic_cost$P_Meh3_cost, na.rm = TRUE), by = 5),
    expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +

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




Rejector_Matrix_Economic_cost_Plot_Map <- egg::ggarrange(
  Clay_economic_cost_plot,
  SOM_economic_cost_plot,
  K_economic_cost_plot,
  P_economic_cost_plot, ncol=2,nrow=2)

ggsave("Figures/Fig4.jpg", plot = Rejector_Matrix_Economic_cost_Plot_Map, width = 3000*1.5, height = 2600*1.5, units = "px", dpi = 300)










