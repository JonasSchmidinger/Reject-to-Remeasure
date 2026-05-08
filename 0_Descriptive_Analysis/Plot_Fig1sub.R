# In this script, we create the different parts for Fig.1, the rest was constructed through Mircrosoft Powerpoint

library(ggplot2)
library(tidyverse)
library(showtext)
library(sysfonts)

font_add(
  family = "Arial",
  regular = "C:/Windows/Fonts/Arial.ttf",
  italic  = "C:/Windows/Fonts/Ariali.ttf"
)
showtext_auto(enable = TRUE)
sysfonts::font_families()


# CDF
delta <- 10
yhat  <- 25
pi_col <- "#7FA8C9"
pi_line <- "#5B8FB3"
x <- seq(0, 50, length.out = 1000)


sd1 <- delta / qnorm(0.995)
sd2 <- delta / qnorm(0.85)



clean_theme <- theme_bw(base_family = "Arial") +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 26, colour = "black"),
    axis.title = element_text(size = 32, colour = "black"),
    axis.ticks = element_blank(),

    axis.line.x = element_line(
      colour = "black", linewidth = 1,
      arrow = arrow(length = unit(0.2, "cm"), type = "closed")
    ),
    axis.line.y = element_line(
      colour = "black", linewidth = 1,
      arrow = arrow(length = unit(0.2, "cm"), type = "closed")
    ),

    text = element_text(family = "Arial")
  )


sym_yhat  <- "\u0177"  # ŷ
sym_delta <- "\u03B4"  # δ
sym_minus <- "\u2212"  # − (proper minus)
sym_geq   <- "\u2265"  # ≥
sym_mid   <- "\u2223"  # ∣ (better conditional bar than |)



make_cdf_plot <- function(mean, sd, x_grid) {
  lower <- mean - delta
  upper <- mean + delta

  df <- data.frame(
    x = x_grid,
    y = pnorm(x_grid, mean, sd)
  )

  df_left  <- df %>% filter(x <= lower)
  df_right <- df %>% filter(x >= upper)

  lbl_yhat      <- sym_yhat
  lbl_yhat_m    <- paste0(sym_yhat, " ", sym_minus, " ", sym_delta)
  lbl_yhat_p    <- paste0(sym_yhat, " + ", sym_delta)


  ggplot(df, aes(x, y)) +

    geom_ribbon(
      data = df_left,
      aes(ymin = 0, ymax = y),
      fill = pi_line, alpha = 0.6
    ) +

    geom_ribbon(
      data = df_right,
      aes(ymin = y, ymax = 1),
      fill = pi_line, alpha = 0.6
    ) +

    geom_line(colour = "black", linewidth = 0.8) +

    geom_segment(
      x = mean, xend = mean,
      y = 0, yend = 1,
      colour = "black",
      linewidth = 1
    ) +
    geom_segment(
      x = lower, xend = lower,
      y = 0, yend = 1,
      linetype = "dotted",
      linewidth = 1.2
    ) +
    geom_segment(
      x = upper, xend = upper,
      y = 0, yend = 1,
      linetype = "dotted",
      linewidth = 1.2
    ) +

    scale_x_continuous(
      limits = c(0, 53),
      breaks = seq(0, 50, 25),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1.02),
      breaks = c(0, 0.5, 1),
      expand = c(0, 0)
    ) +

    clean_theme +

  annotate(
    "text",
    x = 27, y = 0.4,
    label = lbl_yhat,
    family = "Arial",
    fontface = "italic",
    size = 11
  ) +
    annotate(
      "text",
      x = 10, y = 0.4,
      label = lbl_yhat_m,
      family = "Arial",
      fontface = "italic",
      size = 10
    ) +
    annotate(
      "text",
      x = 40, y = 0.4,
      label = lbl_yhat_p,
      family = "Arial",
      fontface = "italic",
      size = 10
    ) +

  labs(
    x = expression(italic(y)),
    y = expression(plain("Pr") * "(" * italic(Y) <= italic(y) * ")")
  )
}


p1_cdf <- make_cdf_plot(yhat, sd1, x)
p2_cdf <- make_cdf_plot(yhat, sd2, x)

p1_cdf
p2_cdf


ggsave("Figures/Fig1_cdf_1.jpg", plot = p1_cdf, width = 800, height = 450, units = "px", dpi = 200)
ggsave("Figures/Fig1_cdf_2.jpg", plot = p2_cdf, width = 800, height = 450, units = "px", dpi = 200)





Full_Dataset <- read.csv("Data/Full_Dataset.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
Full_Dataset

Spectra <- Full_Dataset[-c(2:11)]
Spectra

wl_cols <- grep("^wl_", names(Spectra), value = TRUE)

Spectra_long <- Spectra %>%
  pivot_longer(
    cols = all_of(wl_cols),
    names_to = "wavelength",
    values_to = "reflectance",
    names_prefix = "wl_",                              # remove prefix
    names_transform = list(wavelength = as.numeric)    # now numeric
  )


Spectra_long_shift <- Spectra_long %>%
  mutate(
    reflectance_shift = ifelse(ID == 636,
                               reflectance - 0.08,
                               reflectance)
  )


spectra_plot_uncertain <- ggplot() +

  # --- all spectra EXCEPT 636 in grey ---
  geom_line(
    data = Spectra_long_shift %>% filter(ID != 636),
    aes(x = wavelength, y = reflectance_shift, group = ID),
    colour = "grey60",
    alpha = 0.15,
    linewidth = 0.4
  ) +

  # --- ONLY ID 636 in blue ---
  geom_line(
    data = Spectra_long_shift %>% filter(ID == 636),
    aes(x = wavelength, y = reflectance_shift, group = ID),
    colour = "black",
    linewidth = 1.3
  ) +

  labs(
    x = "Wavelength",
    y = "Reflectance"
  ) +

  theme_bw() +

  theme(
    # --- remove numbers and ticks ---
    axis.text  = element_blank(),
    axis.ticks = element_blank(),

    # --- axis title style ---
    axis.title = element_text(size = 32),

    # --- arrows on axes ---
    axis.line.x = element_line(
      colour = "black",
      linewidth = 1,
      arrow = arrow(length = unit(0.2, "cm"), type = "closed")
    ),

    axis.line.y = element_line(
      colour = "black",
      linewidth = 1,
      arrow = arrow(length = unit(0.2, "cm"), type = "closed")
    ),

    panel.border = element_blank(),
    panel.grid   = element_blank(),

    legend.position = "none",

    text = element_text(family = "Arial")
  )

spectra_plot_uncertain
ggsave("Figures/Fig1_Uncertain_Spectra.jpg", plot = spectra_plot_uncertain, width = 800, height = 450, units = "px", dpi = 200)







spectra_plot_certain <- ggplot() +

  geom_line(
    data = Spectra_long %>% filter(ID != 250),
    aes(x = wavelength, y = reflectance, group = ID),
    colour = "grey60",
    alpha = 0.15,
    linewidth = 0.4
  ) +

  geom_line(
    data = Spectra_long %>% filter(ID == 250),
    aes(x = wavelength, y = reflectance, group = ID),
    colour = "black",
    linewidth = 1.3
  ) +

  labs(
    x = "Wavelength",
    y = "Reflectance"
  ) +

  theme_bw() +

  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),

    axis.title = element_text(size = 32),

    axis.line.x = element_line(
      colour = "black",
      linewidth = 1,
      arrow = arrow(length = unit(0.2, "cm"), type = "closed")
    ),

    axis.line.y = element_line(
      colour = "black",
      linewidth = 1,
      arrow = arrow(length = unit(0.2, "cm"), type = "closed")
    ),

    panel.border = element_blank(),
    panel.grid   = element_blank(),

    legend.position = "none",

    text = element_text(family = "Arial")
  )

spectra_plot_certain
ggsave("Figures/Fig1_Certain_Spectra.jpg", plot = spectra_plot_certain, width = 800, height = 450, units = "px", dpi = 200)
