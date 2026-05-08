# In this script, we create Fig.S1 and Fig.S2

library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(ggpubr)
library(prospectr)
library(distfromq)
library(grid)

library(showtext)
library(sysfonts)
library(ggtext)


font_add(
  family = "Arial",
  regular = "C:/Windows/Fonts/arial.ttf",
  italic  = "C:/Windows/Fonts/ariali.ttf",
  bold    = "C:/Windows/Fonts/arialbd.ttf",
  bolditalic = "C:/Windows/Fonts/arialbi.ttf"
)



quantiles <-
  c(
    0.001,
    round(seq(0.005, 0.995, by = 0.005), 3),
    0.999
  )





SOM_true<- read.csv(file = "Results/TabPFN/Outer/SOM_true.csv")
Clay_true<- read.csv(file = "Results/TabPFN/Outer/Clay_true.csv")
K_Meh3_true<-  read.csv(file = "Results/TabPFN/Outer/K_Meh3_true.csv")
P_Meh3_true<-  read.csv(file = "Results/TabPFN/Outer/P_Meh3_true.csv")
Total_C_true<- read.csv(file = "Results/TabPFN/Outer/Total_C_true.csv")
Sand_true<- read.csv(file = "Results/TabPFN/Outer/Sand_true.csv")
pH_SMP_true<-  read.csv(file = "Results/TabPFN/Outer/pH_SMP_true.csv")
Al_Meh3_true<- read.csv(file = "Results/TabPFN/Outer/Al_Meh3_true.csv")


Total_C_pred<- read.csv(file = "Results/TabPFN/Outer/Total_C_pred.csv")
Sand_pred<- read.csv(file = "Results/TabICL/Outer/Sand_pred.csv")
pH_SMP_pred<-  read.csv(file = "Results/TabPFN/Outer/pH_SMP_pred.csv")
Al_Meh3_pred<- read.csv(file = "Results/TabICL/Outer/Al_Meh3_pred.csv")
SOM_pred<- read.csv(file = "Results/TabPFN/Outer/SOM_pred.csv")
Clay_pred<- read.csv(file = "Results/TabICL/Outer/Clay_pred.csv")
K_Meh3_pred<-  read.csv(file = "Results/TabPFN/Outer/K_Meh3_pred.csv")
P_Meh3_pred<-  read.csv(file = "Results/TabPFN/Outer/P_Meh3_pred.csv")


Total_C_pred_quantiles<- read.csv(file = "Results/TabPFN/Outer/Total_C_pred_quantiles.csv")
Sand_pred_quantiles<- read.csv(file = "Results/TabICL/Outer/Sand_pred_quantiles.csv")
pH_SMP_pred_quantiles<-  read.csv(file = "Results/TabPFN/Outer/pH_SMP_pred_quantiles.csv")
Al_Meh3_pred_quantiles<- read.csv(file = "Results/TabICL/Outer/Al_Meh3_pred_quantiles.csv")
SOM_pred_quantiles<- read.csv(file = "Results/TabPFN/Outer/SOM_pred_quantiles.csv")
Clay_pred_quantiles<- read.csv(file = "Results/TabICL/Outer/Clay_pred_quantiles.csv")
K_Meh3_pred_quantiles<-  read.csv(file = "Results/TabPFN/Outer/K_Meh3_pred_quantiles.csv")
P_Meh3_pred_quantiles<-  read.csv(file = "Results/TabPFN/Outer/P_Meh3_pred_quantiles.csv")









quantile_coverage_function <- function(test, quantile_predictions, quantiles) {
  if (ncol(quantile_predictions) != length(quantiles)) {
    stop("Number of columns in quantile_predictions must equal length(quantiles).")
  }

  n <- length(test)

  coverage <- sapply(seq_along(quantiles), function(j) {
    mean(test <= quantile_predictions[, j]) * 100
  })

  in_counts <- sapply(seq_along(quantiles), function(j) {
    sum(test <= quantile_predictions[, j])
  })

  out_counts <- n - in_counts

  coverage.df <- data.frame(
    quantile = quantiles,
    coverage = coverage,
    within = in_counts,
    out = out_counts
  )

  return(coverage.df)
}


# --- Only keep relevant mapping ---
true_list <- list(
  SOM     = SOM_true[[1]],
  Clay    = Clay_true[[1]],
  K_Meh3  = K_Meh3_true[[1]],
  P_Meh3  = P_Meh3_true[[1]],
  Total_C = Total_C_true[[1]],
  Sand    = Sand_true[[1]],
  Al_Meh3 = Al_Meh3_true[[1]],
  pH_SMP  = pH_SMP_true[[1]]
)

# --- Only the models you actually use ---
quantile_pred_list <- list(
  # TabPFN
  SOM     = SOM_pred_quantiles,
  K_Meh3  = K_Meh3_pred_quantiles,
  P_Meh3  = P_Meh3_pred_quantiles,
  Total_C = Total_C_pred_quantiles,
  pH_SMP  = pH_SMP_pred_quantiles,

  # TabICL
  Clay    = Clay_pred_quantiles,
  Sand    = Sand_pred_quantiles,
  Al_Meh3 = Al_Meh3_pred_quantiles
)

model_map <- c(
  SOM = "TabPFN",
  K_Meh3 = "TabPFN",
  P_Meh3 = "TabPFN",
  Total_C = "TabPFN",
  pH_SMP = "TabPFN",
  Clay = "TabICL",
  Sand = "TabICL",
  Al_Meh3 = "TabICL"
)

# --- Compute ONLY subset ---
coverage_all <- lapply(names(quantile_pred_list), function(nm) {

  res <- quantile_coverage_function(
    test = true_list[[nm]],
    quantile_predictions = quantile_pred_list[[nm]],
    quantiles = quantiles
  )

  res$Property <- nm
  res$Model <- model_map[[nm]]

  res
}) %>% bind_rows()




# --- Clean naming ---
coverage_all <- coverage_all %>%
  mutate(
    Property = recode(Property,
                      "K_Meh3" = "K",
                      "P_Meh3" = "P",
                      "Total_C" = "TC",
                      "Al_Meh3" = "Al",
                      "pH_SMP" = "pH"),
    Property = factor(Property,
                      levels = c("Clay","SOM","K","P","Sand","TC","Al","pH"))
  )
coverage_all <- coverage_all %>%
  mutate(
    facet_label = paste0(Property, " (", Model, ")")
  )




QCP_Plots <- ggplot(
  coverage_all,
  aes(x = quantile * 100,
      y = coverage)
) +

  geom_line(
    linewidth = 2,
    alpha = 0.8,
    colour = "#1F4E79"
  ) +

  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 1.2,
    colour = "black"
  ) +

  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.02)),
    limits = c(0, 100)
  ) +

  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.02)),
    limits = c(0, 100)
  ) +

  labs(
    x = "Nominal Quantile Level (%)",
    y = "QCP (%)"
  ) +

  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 20, colour = "black"),
    legend.position = "none",
    panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
    axis.line = element_blank(),
    strip.text = element_text(size = 18, colour = "black", face = "bold"),
    strip.background = element_rect(fill = "grey95", linewidth = 0.5)
  ) +

  facet_wrap(~ facet_label, ncol = 2, nrow = 4)


QCP_Plots

ggsave("Figures/FigS1.jpg", plot = QCP_Plots,  width = 3000*1.5,
       height = 4500*1.5, units = "px", dpi = 300)





###









selected_quantiles <- c(0.025, 0.975)
n_bins <- 5


calibration_deviation_function <- function(y_true,
                                           quantile_preds,
                                           quantiles,
                                           selected_q = c(0.025, 0.975),
                                           n_bins = 5,
                                           round_digits = 0){

  breaks <- quantile(y_true, probs = seq(0,1,length.out = n_bins+1), na.rm = TRUE)

  bin <- cut(y_true, breaks = breaks, include.lowest = TRUE, labels = FALSE)

  bin_labels <- paste0(
    round(breaks[-length(breaks)], round_digits),
    "–",
    round(breaks[-1], round_digits)
  )

  results <- list()

  for(q in selected_q){

    q_index <- which.min(abs(quantiles - q))
    pred_q <- quantile_preds[, q_index]

    df <- data.frame(y = y_true, pred = pred_q, bin = bin)

    res <- df %>%
      dplyr::group_by(bin) %>%
      dplyr::summarise(coverage = mean(y <= pred), .groups = "drop")

    res$deviation <- (res$coverage - q) * 100
    res$quantile <- q
    res$bin_label <- factor(bin_labels[res$bin], levels = bin_labels)

    results[[as.character(q)]] <- res
  }

  dplyr::bind_rows(results)
}


theme_calibration <- theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 20, colour = "black"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
    strip.text = element_text(size = 18, face = "bold"),
    strip.background = element_rect(fill = "grey95")
  )


panel_config <- data.frame(
  prop  = c("Clay","SOM","K","P","Sand","TC","Al","pH"),
  model = c("TabICL","TabPFN","TabPFN","TabPFN","TabICL","TabPFN","TabICL","TabPFN"),
  label = c(
    "bold(Clay~(g~kg^{-1})*' - TabICL')",
    "bold(SOM~(g~kg^{-1})*' - TabPFN')",
    "bold(K~(mg~kg^{-1})*' - TabPFN')",
    "bold(P~(mg~kg^{-1})*' - TabPFN')",
    "bold(Sand~(g~kg^{-1})*' - TabICL')",
    "bold(TC~(g~kg^{-1})*' - TabPFN')",
    "bold(Al~(mg~kg^{-1})*' - TabICL')",
    "bold(pH~' - TabPFN')"
  )
)


plot_list <- lapply(seq_len(nrow(panel_config)), function(i){

  prop_name  <- panel_config$prop[i]
  model_name <- panel_config$model[i]
  facet_lab  <- panel_config$label[i]

  is_left <- i %% 2 == 1

  prop_raw <- dplyr::recode(prop_name,
                            "K"="K_Meh3",
                            "P"="P_Meh3",
                            "TC"="Total_C",
                            "Al"="Al_Meh3",
                            "pH"="pH_SMP")

  y_true <- true_list[[prop_raw]]
  pred_df <- quantile_pred_list[[prop_raw]]

  if(prop_name %in% c("Clay","SOM","Sand","TC")){
    y_true  <- y_true * 10
    pred_df <- pred_df * 10
  }


  round_digits <- if(prop_name == "pH") 1 else 0

  df <- calibration_deviation_function(
    y_true,
    pred_df,
    quantiles,
    selected_q = selected_quantiles,
    n_bins = n_bins,
    round_digits = round_digits
  )

  df$quantile <- factor(df$quantile,
                        levels = c(0.025,0.975),
                        labels = c("2.5","97.5"))

  df$facet_label <- facet_lab

  y_title <- if(prop_name == "K") "Coverage Error (%)" else NULL
  x_title <- if(prop_name == "Al") "Target Range (Y)" else NULL



  p <- ggplot(df,
              aes(x = bin,
                  y = deviation,
                  colour = quantile,
                  group = quantile)) +

    geom_segment(
      aes(x = bin - 0.3,
          xend = bin + 0.3,
          yend = deviation),
      linewidth = 2,
      alpha = 0.7
    ) +

    geom_point(size = 4, alpha = 0.7) +

    geom_hline(yintercept = 0,
               colour = "black",
               linewidth = 1) +

    scale_colour_manual(
      values = c("purple3","#1B9E77"),
      name = "Nominal Quantile Level (%):"
    ) +

    scale_y_continuous(limits = c(-20, 20)) +

    scale_x_continuous(
      breaks = sort(unique(df$bin)),
      labels = levels(df$bin_label),
      expand = expansion(mult = c(0.05, 0.05))
    ) +

    labs(x = x_title, y = y_title) +

    theme_calibration +

    facet_wrap(~facet_label, labeller = label_parsed)



  if(!is_left){
    p <- p + theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank()
    )
  }



  if(prop_name == "K"){
    p <- p + theme(axis.title.y = element_text(hjust = -0.5))
  }

  if(prop_name == "Al"){
    p <- p + theme(axis.title.x = element_text(hjust = 1.28))
  }



  if(prop_name == "SOM"){
    p <- p + theme(
      legend.position = c(-0.45, 1.25),
      legend.justification = c(0, 1),
      legend.direction = "horizontal",
      legend.background = element_rect(fill = "transparent"),
      plot.margin = margin(t = 60, r = 5, b = 5, l = 5)
    ) +
      guides(
        colour = guide_legend(
          nrow = 1,
          byrow = TRUE,
          override.aes = list(
            size = 6,
            linewidth = 3
          ),
          keywidth = unit(2, "cm"),
          keyheight = unit(0.8, "cm")
        )
      )
  } else {
    p <- p + theme(legend.position = "none")
  }

  return(p)
})


Calibration_Final <- do.call(
  egg::ggarrange,
  c(plot_list, list(ncol = 2, nrow = 4))
)

Calibration_Final

ggsave(
  "Figures/FigS2.jpg",
  plot = Calibration_Final,
  width = 3000 * (300/200),
  height = 4600 * (300/200),
  units = "px",
  dpi = 300
)

