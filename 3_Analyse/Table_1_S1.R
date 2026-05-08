library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(ggpubr)
library(prospectr)
library(ggpointdensity)


R2 <- function(y_hat, y) {
  1 - sum((y - y_hat)^2) / sum((y - mean(y))^2)
}

MAE <- function(y_hat, y) {
  mean(abs((y - y_hat)))
}

quantiles <- c(
  0.001,
  round(seq(0.005, 0.995, by = 0.005), 3),
  0.999
)

quantile_coverage_function <- function(test, quantile_predictions, quantiles) {
  coverage <- sapply(seq_along(quantiles), function(j) {
    mean(test <= quantile_predictions[, j]) * 100
  })

  data.frame(
    quantile = quantiles,
    coverage = coverage
  )
}

crps_step_matrix <- function(y, pred_mat, quantiles) {
  w <- diff(c(0, quantiles, 1))

  sapply(seq_len(nrow(pred_mat)), function(i) {
    q_row <- as.numeric(pred_mat[i, ])
    q_ext <- c(q_row, tail(q_row, 1))

    term1 <- sum(w * abs(q_ext - y[i]))
    pairwise_mat <- abs(outer(q_ext, q_ext, "-"))
    term2 <- 0.5 * sum(w * (pairwise_mat %*% w))

    term1 - term2
  })
}

compute_mpiw_95 <- function(quant_df) {
  mean(quant_df$q0.975 - quant_df$q0.025, na.rm = TRUE)
}


compute_metrics <- function(true_df, pred_df, quant_df, quantiles, target_name) {

  y     <- true_df[[1]]
  y_hat <- pred_df[[1]]
  qmat  <- as.matrix(quant_df)

  if (target_name %in% c("Clay", "SOM", "Sand", "Total_C")) {
    y     <- y * 10
    y_hat <- y_hat * 10
    qmat  <- qmat * 10
    quant_df <- as.data.frame(qmat)
  }


  r2_val  <- R2(y_hat, y)
  mae_val <- MAE(y_hat, y)

  coverage_df <- quantile_coverage_function(
    test = y,
    quantile_predictions = as.matrix(quant_df[, c("q0.025", "q0.975")]),
    quantiles = c(0.025, 0.975)
  )

  QCP0.025 <- coverage_df$coverage[coverage_df$quantile == 0.025]
  QCP0.975 <- coverage_df$coverage[coverage_df$quantile == 0.975]

  MPIW_95 <- mean(quant_df$q0.975 - quant_df$q0.025, na.rm = TRUE)

  CRPS <- mean(crps_step_matrix(y, qmat, quantiles), na.rm = TRUE)

  data.frame(
    Target = NA,
    Model = NA,
    MPIW_95 = MPIW_95,
    QCP0.025 = QCP0.025,
    QCP0.975 = QCP0.975,
    CRPS = CRPS,
    R2 = r2_val,
    Mean_AE = mae_val
  )
}


Clay_true     <- read.csv("Results/TabPFN/Outer/Clay_true.csv")
SOM_true      <- read.csv("Results/TabPFN/Outer/SOM_true.csv")
K_Meh3_true   <- read.csv("Results/TabPFN/Outer/K_Meh3_true.csv")
P_Meh3_true   <- read.csv("Results/TabPFN/Outer/P_Meh3_true.csv")
Sand_true     <- read.csv("Results/TabPFN/Outer/Sand_true.csv")
Total_C_true  <- read.csv("Results/TabPFN/Outer/Total_C_true.csv")
Al_Meh3_true  <- read.csv("Results/TabPFN/Outer/Al_Meh3_true.csv")
pH_SMP_true   <- read.csv("Results/TabPFN/Outer/pH_SMP_true.csv")


load_model <- function(model) {
  list(
    Clay_pred = read.csv(paste0("Results/", model, "/Outer/Clay_pred.csv")),
    Clay_q    = read.csv(paste0("Results/", model, "/Outer/Clay_pred_quantiles.csv")),

    SOM_pred = read.csv(paste0("Results/", model, "/Outer/SOM_pred.csv")),
    SOM_q    = read.csv(paste0("Results/", model, "/Outer/SOM_pred_quantiles.csv")),

    K_pred = read.csv(paste0("Results/", model, "/Outer/K_Meh3_pred.csv")),
    K_q    = read.csv(paste0("Results/", model, "/Outer/K_Meh3_pred_quantiles.csv")),

    P_pred = read.csv(paste0("Results/", model, "/Outer/P_Meh3_pred.csv")),
    P_q    = read.csv(paste0("Results/", model, "/Outer/P_Meh3_pred_quantiles.csv")),

    Sand_pred = read.csv(paste0("Results/", model, "/Outer/Sand_pred.csv")),
    Sand_q    = read.csv(paste0("Results/", model, "/Outer/Sand_pred_quantiles.csv")),

    TC_pred = read.csv(paste0("Results/", model, "/Outer/Total_C_pred.csv")),
    TC_q    = read.csv(paste0("Results/", model, "/Outer/Total_C_pred_quantiles.csv")),

    Al_pred = read.csv(paste0("Results/", model, "/Outer/Al_Meh3_pred.csv")),
    Al_q    = read.csv(paste0("Results/", model, "/Outer/Al_Meh3_pred_quantiles.csv")),

    pH_pred = read.csv(paste0("Results/", model, "/Outer/pH_SMP_pred.csv")),
    pH_q    = read.csv(paste0("Results/", model, "/Outer/pH_SMP_pred_quantiles.csv"))
  )
}

models <- list(
  PLSQR = load_model("PLSQR"),
  QRF   = load_model("QRF"),
  TabPFN = load_model("TabPFN"),
  TabICL = load_model("TabICL")
)


results_list <- list()

for (model_name in names(models)) {

  m <- models[[model_name]]

  add_row <- function(target, true, pred, quant) {
    df <- compute_metrics(true, pred, quant, quantiles, target_name = target)
    df$Target <- target
    df$Model  <- model_name
    df
  }

  results_list[[length(results_list)+1]] <- add_row("Clay",     Clay_true,   m$Clay_pred, m$Clay_q)
  results_list[[length(results_list)+1]] <- add_row("SOM",      SOM_true,    m$SOM_pred,  m$SOM_q)
  results_list[[length(results_list)+1]] <- add_row("K_Meh3",   K_Meh3_true, m$K_pred,    m$K_q)
  results_list[[length(results_list)+1]] <- add_row("P_Meh3",   P_Meh3_true, m$P_pred,    m$P_q)
  results_list[[length(results_list)+1]] <- add_row("Sand",     Sand_true,   m$Sand_pred, m$Sand_q)
  results_list[[length(results_list)+1]] <- add_row("Total_C",  Total_C_true,m$TC_pred,   m$TC_q)
  results_list[[length(results_list)+1]] <- add_row("Al_Meh3",  Al_Meh3_true,m$Al_pred,   m$Al_q)
  results_list[[length(results_list)+1]] <- add_row("pH_SMP",   pH_SMP_true, m$pH_pred,   m$pH_q)
}

results_table <- do.call(rbind, results_list)

results_table <- results_table[, c(
  "Target", "Model", "MPIW_95", "QCP0.025", "QCP0.975", "CRPS", "R2", "Mean_AE"
)]

results_table <- results_table %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

target_order <- c(
  "Clay", "SOM", "K_Meh3", "P_Meh3",
  "Sand", "Total_C", "Al_Meh3", "pH_SMP"
)

model_order <- c("PLSQR", "QRF", "TabPFN", "TabICL")

results_table <- results_table %>%
  mutate(
    Target = factor(Target, levels = target_order),
    Model  = factor(Model,  levels = model_order)
  ) %>%
  arrange(Target, Model) %>%
  mutate(
    Target = as.character(Target),
    Model  = as.character(Model)
  )

# save
#write.csv(results_table, "Results/results_table.csv", row.names = FALSE)

results_table
