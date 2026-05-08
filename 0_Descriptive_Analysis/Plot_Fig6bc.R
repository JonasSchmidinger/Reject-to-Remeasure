# In this script, we create Fig. 6b and Fig. 6c. Since we do not have the rights to share coordinates, we cannot share information about Fig. 6a

library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(showtext)
library(sysfonts)

font_add(
  family = "Arial",
  regular = "C:/Windows/Fonts/Arial.ttf",
  italic  = "C:/Windows/Fonts/Ariali.ttf"
)
showtext_auto(enable = TRUE)
sysfonts::font_families()



Dataset <- read.csv("Data/Full_Dataset.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
Target_Properties <-Dataset[c("Clay","SOM","K_Meh3","P_Meh3","Sand","Total_C","Al_Meh3","pH_SMP")]

Target_Properties <- Dataset %>%
  select(Clay, SOM, K_Meh3, P_Meh3, Sand, Total_C, Al_Meh3, pH_SMP) %>%
  rename(
    K  = K_Meh3,
    P  = P_Meh3,
    TC = Total_C,
    Al = Al_Meh3,
    pH = pH_SMP
  )



cor_mat <- cor(Target_Properties, use = "pairwise.complete.obs")

cor_df <- as.data.frame(as.table(cor_mat)) %>%
  rename(Var1 = Var1, Var2 = Var2, Correlation = Freq)

cor_df <- cor_df %>%
  mutate(
    Var1 = factor(Var1, levels = colnames(cor_mat)),
    Var2 = factor(Var2, levels = colnames(cor_mat))
  ) %>%
  filter(as.numeric(Var1) > as.numeric(Var2)) %>%
  droplevels()   # ← THIS removes the empty "Clay" row

corrplot <- ggplot(cor_df, aes(x = Var2, y = Var1, fill = Correlation)) +

  geom_tile(color = "white") +

  geom_text(
    aes(label = sprintf("%.2f", Correlation)),
    color = "black",
    size = 6
  ) +

  scale_fill_gradient2(
    low = "#6A00A8",
    mid = "white",
    high = "#00A86B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson's Correlation",
    guide = guide_colourbar(
      direction = "horizontal",
      title.position = "bottom",
      title.hjust = 0.5,
      ticks = TRUE,
      ticks.colour = "black",
      frame.colour = "black",
      frame.linewidth = 0.3
    )
  ) +

  coord_fixed() +

  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(
    limits = rev(levels(cor_df$Var1)),
    expand = c(0, 0)
  ) +

  theme_bw(base_family = "Arial") +
  theme(
    text = element_text(family = "Arial"),

    panel.border = element_blank(),
    panel.grid = element_blank(),

    axis.title = element_blank(),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 28,
      color= "black"
    ),
    axis.text.y = element_text(size = 28,color= "black"),

    axis.ticks = element_line(color = "black"),

    legend.position = "bottom",
    legend.direction = "horizontal",

    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16),

    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(0.4, "cm"),

    legend.background = element_blank()
  )

corrplot

ggsave("Figures/Fig6b.jpg", plot = corrplot, width = 900, height = 1200, units = "px", dpi = 200)


#



Spectra <- Dataset[-c(2:20)]



wl_cols <- grep("^wl_", names(Spectra), value = TRUE)

Spectra_long <- Spectra %>%
  pivot_longer(
    cols = all_of(wl_cols),
    names_to = "wavelength",
    values_to = "reflectance",
    names_prefix = "wl_",
    names_transform = list(wavelength = as.numeric)
  )

Spectra_long


spectra_plot <- ggplot() +

  # --- all spectra EXCEPT 636 in grey ---
  geom_line(
    data = Spectra_long,
    aes(x = wavelength , y = reflectance*100, group = ID),
    colour = "grey40",
    alpha = 0.2,
    linewidth = 0.4
  ) +

  labs(
    x = "Wavelength (nm)",
    y = "Reflectance (%)"
  ) +

  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),limits = c(-0.5,80.5)) +

  theme_bw() +

  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 24,color= "black"),
    axis.title = element_text(size = 28,color= "black"),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 12.5, b = 5, l = 5))


spectra_plot
ggsave("Figures/Fig6c.jpg", plot = spectra_plot, width = 900, height = 1200, units = "px", dpi = 200)



