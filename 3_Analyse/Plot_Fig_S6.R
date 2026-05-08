# In this script, we create Fig.S6

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


###

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

###

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
# Sand

Sand_abs_buffer      <- 120
Sand_rel_factor      <- 0
Sand_threshold_alpha <- 0.05



Sand_p_value <- compute_p_exceed_hybrid(
  quantile_matrix = as.matrix(Sand_pred_quantiles) * 10,
  pred_values     = Sand_pred$Sand_pred * 10,
  a               = Sand_abs_buffer,
  b               = Sand_rel_factor,
  quantile_probs  = quantiles
)



Sand_error_p.df <- data.frame(
  pred    = as.numeric(Sand_pred$Sand_pred) * 10,
  true    = as.numeric(Sand_true$Sand_target) * 10,
  p_value = Sand_p_value$p_value
) %>%
  mutate(
    error = abs(pred - true),

    delta_true = Sand_abs_buffer + Sand_rel_factor * true,

    is_wrong = error > delta_true,

    is_rejected = p_value > Sand_threshold_alpha
  )



Sand_acceptance_rate <- mean(!Sand_error_p.df$is_rejected)




visnir_costs  <- seq(0.75, 10.7, by = 0.025)
texture_costs <- seq(33, 56.5, by = 0.025) # Slightly adjusted range for better visibility (done for all ranges !)

Sand_economic_cost <- expand.grid(
  visnir_cost  = visnir_costs,
  texture_cost = texture_costs
)


Sand_economic_cost <- Sand_economic_cost %>%
  mutate(

    baseline_cost = texture_cost,

    strategy_cost = visnir_cost +
      (1 - Sand_acceptance_rate) * texture_cost,

    cost_difference = strategy_cost - baseline_cost
  )


Sand_acceptance_percent <- round(Sand_acceptance_rate * 100, 1)

Sand_economic_cost$title <- paste0(
  "bold(Sand~'-'~'",
  Sand_acceptance_percent,
  "% Acceptance Rate')"
)




Sand_economic_cost <- Sand_economic_cost %>%
  mutate(
    visnir_cost  = as.numeric(visnir_cost),
    texture_cost = as.numeric(texture_cost),
    cost_difference  = as.numeric(cost_difference)
  )


contour_levels <- seq(-40, 40, by = 5)


Sand_economic_cost_plot <- ggplot(
  Sand_economic_cost,
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

  # geom_contour(
  #   data = Sand_economic_cost,
  #   aes(
  #     x = visnir_cost,
  #     y = texture_cost,
  #     z = cost_difference,
  #     group = 1
  #   ),
  #   breaks = 0,
  #   colour = "black",
  #   linewidth = 0.8,
  #   inherit.aes = FALSE
  # ) +

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
    breaks = seq(35, max(Sand_economic_cost$texture_cost, na.rm = TRUE), by = 10),
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

Sand_economic_cost_plot

#----------------------------------------------------------------------------------------------------------------
# Total_C


Total_C_abs_buffer      <- 10
Total_C_rel_factor      <- 0
Total_C_threshold_alpha <- 0.05


Total_C_p_value <- compute_p_exceed_hybrid(
  quantile_matrix = as.matrix(Total_C_pred_quantiles) * 10,
  pred_values     = Total_C_pred$Total_C_pred * 10,
  a               = Total_C_abs_buffer,
  b               = Total_C_rel_factor,
  quantile_probs  = quantiles
)


Total_C_error_p.df <- data.frame(
  pred    = as.numeric(Total_C_pred$Total_C_pred) * 10,
  true    = as.numeric(Total_C_true$Total_C_target) * 10,
  p_value = Total_C_p_value$p_value
) %>%
  mutate(
    error = abs(pred - true),

    delta_true = Total_C_abs_buffer + Total_C_rel_factor * true,

    is_wrong = error > delta_true,

    is_rejected = p_value > Total_C_threshold_alpha
  )



Total_C_acceptance_rate <- mean(!Total_C_error_p.df$is_rejected)



visnir_costs  <- seq(0.75, 10.7, by = 0.025)
Total_C_costs <- seq(14.5, 26, by = 0.025)

Total_C_economic_cost <- expand.grid(
  visnir_cost  = visnir_costs,
  Total_C_cost = Total_C_costs
)


Total_C_economic_cost <- Total_C_economic_cost %>%
  mutate(

    baseline_cost = Total_C_cost,

    strategy_cost = visnir_cost +
      (1 - Total_C_acceptance_rate) * Total_C_cost,

    cost_difference = strategy_cost - baseline_cost
  )


Total_C_acceptance_percent <- round(Total_C_acceptance_rate * 100, 1)

Total_C_economic_cost$title <- paste0(
  "bold(TC~'-'~'",
  Total_C_acceptance_percent,
  "% Acceptance Rate')"
)



Total_C_economic_cost <- Total_C_economic_cost %>%
  mutate(
    visnir_cost  = as.numeric(visnir_cost),
    Total_C_cost = as.numeric(Total_C_cost),
    cost_difference  = as.numeric(cost_difference)
  )


Total_C_economic_cost_plot <-ggplot(
  Total_C_economic_cost,
  aes(
    x = visnir_cost,
    y = Total_C_cost,
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
    breaks = seq(15, max(Total_C_economic_cost$Total_C_cost, na.rm = TRUE), by = 5),
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




#----------------------------------------------------------------------------------------------------------------
# Al



Al_Meh3_abs_buffer      <- 250
Al_Meh3_rel_factor      <- 0
Al_Meh3_threshold_alpha <- 0.05


Al_Meh3_p_value <- compute_p_exceed_hybrid(
  quantile_matrix = as.matrix(Al_Meh3_pred_quantiles),
  pred_values     = Al_Meh3_pred$Al_Meh3_pred,
  a               = Al_Meh3_abs_buffer,
  b               = Al_Meh3_rel_factor,
  quantile_probs  = quantiles
)


Al_Meh3_error_p.df <- data.frame(
  pred    = as.numeric(Al_Meh3_pred$Al_Meh3_pred),
  true    = as.numeric(Al_Meh3_true$Al_Meh3_target),
  p_value = Al_Meh3_p_value$p_value
) %>%
  mutate(
    error = abs(pred - true),

    delta_true = Al_Meh3_abs_buffer + Al_Meh3_rel_factor * true,

    is_wrong = error > delta_true,

    is_rejected = p_value > Al_Meh3_threshold_alpha
  )



Al_Meh3_acceptance_rate <- mean(!Al_Meh3_error_p.df$is_rejected)



visnir_costs  <- seq(0.75, 10.7, by = 0.025)
Al_Meh3_costs <-  seq(10, 22, by = 0.025)

Al_Meh3_economic_cost <- expand.grid(
  visnir_cost  = visnir_costs,
  Al_Meh3_cost = Al_Meh3_costs
)


Al_Meh3_economic_cost <- Al_Meh3_economic_cost %>%
  mutate(

    baseline_cost = Al_Meh3_cost,

    strategy_cost = visnir_cost +
      (1 - Al_Meh3_acceptance_rate) * Al_Meh3_cost,

    cost_difference = strategy_cost - baseline_cost
  )


Al_Meh3_acceptance_percent <- round(Al_Meh3_acceptance_rate * 100, 1)

Al_Meh3_economic_cost$title <- paste0(
  "bold(Al~'-'~'",
  Al_Meh3_acceptance_percent,
  "% Acceptance Rate')"
)



Al_Meh3_economic_cost <- Al_Meh3_economic_cost %>%
  mutate(
    visnir_cost  = as.numeric(visnir_cost),
    Al_Meh3_cost = as.numeric(Al_Meh3_cost),
    cost_difference  = as.numeric(cost_difference)
  )


K_economic_cost_plot <-ggplot(
  Al_Meh3_economic_cost,
  aes(
    x = visnir_cost,
    y = Al_Meh3_cost,
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
    name = expression("Cost Saving with VisNIRS" ~ "(" * "$" * "CAD Sample"^{-1} * ")"))+

  scale_y_continuous(expand = c(0,0),
                     breaks = seq(10, max(Al_Meh3_economic_cost$Al_Meh3_cost, na.rm = TRUE), by = 5)) +
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






#----------------------------------------------------------------------------------------------------------------
# pH




pH_SMP_abs_buffer      <- 0.25
pH_SMP_rel_factor      <- 0
pH_SMP_threshold_alpha <- 0.05



pH_SMP_p_value <- compute_p_exceed_hybrid(
  quantile_matrix = as.matrix(pH_SMP_pred_quantiles),
  pred_values     = pH_SMP_pred$pH_SMP_pred,
  a               = pH_SMP_abs_buffer,
  b               = pH_SMP_rel_factor,
  quantile_probs  = quantiles
)


pH_SMP_error_p.df <- data.frame(
  pred    = as.numeric(pH_SMP_pred$pH_SMP_pred),
  true    = as.numeric(pH_SMP_true$pH_SMP_target),
  p_value = pH_SMP_p_value$p_value
) %>%
  mutate(
    error = abs(pred - true),

    delta_true = pH_SMP_abs_buffer + pH_SMP_rel_factor * true,

    is_wrong = error > delta_true,

    is_rejected = p_value > pH_SMP_threshold_alpha
  )



pH_SMP_acceptance_rate <- mean(!pH_SMP_error_p.df$is_rejected)



visnir_costs  <- seq(0.75, 10.7, by = 0.025)
pH_SMP_costs <- seq(10, 22, by = 0.025)

pH_SMP_economic_cost <- expand.grid(
  visnir_cost  = visnir_costs,
  pH_SMP_cost = pH_SMP_costs
)


pH_SMP_economic_cost <- pH_SMP_economic_cost %>%
  mutate(

    baseline_cost = pH_SMP_cost,

    strategy_cost = visnir_cost +
      (1 - pH_SMP_acceptance_rate) * pH_SMP_cost,

    cost_difference = strategy_cost - baseline_cost
  )


pH_SMP_acceptance_percent <- round(pH_SMP_acceptance_rate * 100, 1)

pH_SMP_economic_cost$title <- paste0(
  "bold(pH~'-'~'",
  pH_SMP_acceptance_percent,
  "% Acceptance Rate')"
)


pH_SMP_economic_cost <- pH_SMP_economic_cost %>%
  mutate(
    visnir_cost  = as.numeric(visnir_cost),
    pH_SMP_cost = as.numeric(pH_SMP_cost),
    cost_difference  = as.numeric(cost_difference)
  )


P_economic_cost_plot <-ggplot(
  pH_SMP_economic_cost,
  aes(
    x = visnir_cost,
    y = pH_SMP_cost,
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
    name = expression("Cost saving with VisNIRS" ~ "(" * "$" * "CAD Measurment"^{-1} * ")"))+

  scale_y_continuous(
    breaks = seq(10, max(pH_SMP_economic_cost$pH_SMP_cost, na.rm = TRUE), by = 5),
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




Rejector_Matrix_Economic_cost_Plot_Map_Appendix <- egg::ggarrange(
  Sand_economic_cost_plot,
  Total_C_economic_cost_plot,
  K_economic_cost_plot,
  P_economic_cost_plot, ncol=2,nrow=2)

ggsave("Figures/FigS6.jpg", plot = Rejector_Matrix_Economic_cost_Plot_Map_Appendix, width = 3000*1.5, height = 2600*1.5, units = "px", dpi = 300)










