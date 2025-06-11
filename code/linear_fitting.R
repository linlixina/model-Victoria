#---------------
# LINEAR FITTING
#---------------

library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

df <- read_excel("data/Age_LGA.xlsx")

columns <- c(
  "Rate(Age<12)",
  "Rate(12<=Age<16)",
  "Rate(16=<V<50)",
  "Rate(16=<UV<50)",
  "Rate(V>=50)",
  "Rate(UV>=50)"
)

df_long <- df %>%
  select(overall_Vacc_Level_average, all_of(columns)) %>%
  pivot_longer(
    cols      = all_of(columns),
    names_to  = "group",
    values_to = "PCR"
  ) %>%
  mutate(PCR = ifelse(PCR == 0, 1e-7, PCR))

y_labels <- list(
  "Rate(Age<12)"     = expression(bold(PI['<12y'])),
  "Rate(12<=Age<16)" = expression(bold(PI['12-15y'])),
  "Rate(16=<V<50)"   = expression(bold(PI['16-49y, V'])),
  "Rate(16=<UV<50)"  = expression(bold(PI['16-49y, U'])),
  "Rate(V>=50)"      = expression(bold(PI['50+y, V'])),
  "Rate(UV>=50)"     = expression(bold(PI['50+y, U']))
)

y_axis_settings <- list(
  "Rate(Age<12)"     = list(breaks = c(-0.1,-0.05,0, 0.05, 0.1), limits = c(-0.1, 0.1)),
  "Rate(12<=Age<16)" = list(breaks = c(-0.1,-0.05,0, 0.05, 0.1), limits = c(-0.1, 0.1)),
  "Rate(16=<V<50)"   = list(breaks = c(-0.02,-0.01,0, 0.01, 0.02), limits = c(-0.02, 0.02)),
  "Rate(16=<UV<50)"  = list(breaks = c(-0.1,-0.05,0, 0.05, 0.1,0.15), limits = c(-0.1, 0.15)),
  "Rate(V>=50)"      = list(breaks = c(-0.02,-0.01,0, 0.01, 0.02), limits = c(-0.02, 0.02)),
  "Rate(UV>=50)"     = list(breaks = c(-0.1,-0.05,0, 0.05, 0.1), limits = c(-0.1, 0.1))
)

plot_list <- lapply(seq_along(columns), function(i) {
  col <- columns[i]
  df_sub <- df_long %>% filter(group == col)
  lm_mod <- lm(PCR ~ overall_Vacc_Level_average, data = df_sub)
  coefs  <- coef(lm_mod)
  intercept <- coefs[1]; slope <- coefs[2]
  r2 <- summary(lm_mod)$r.squared
  eq_label <- sprintf("y = %.4f + %.4f x\nR² = %.4f",
                      intercept, slope, r2)
  
  ggplot(df_sub, aes(x = overall_Vacc_Level_average, y = PCR)) +
    geom_smooth(aes(color = "Linear fit", linetype = "Linear fit"),
                method = "lm", se = FALSE, fullrange = TRUE, size = 1) +
    geom_point(aes(color = "Observed data"), size = 1.5, alpha = 0.7) +
    annotate("text",
             x = 0.4,
             y = y_axis_settings[[col]]$limits[2] * 0.8,
             label = eq_label,
             hjust = 0, size = 3.5) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(
      limits = y_axis_settings[[col]]$limits,
      breaks = y_axis_settings[[col]]$breaks
    ) +
    scale_color_manual(name = NULL,
                       values = c(
                         "Linear fit"    = "#33a02c",
                         "Observed data" = "#7F3B08"
                       )) +
    scale_linetype_manual(name = NULL,
                          values = c("Linear fit" = "solid"),guide = "none") +
    labs(
      title = LETTERS[i],
      x     = "Overall Vaccination Level v",
      y     = y_labels[[col]]
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title         = element_text(size = 14, face = "bold", hjust = -0.1),
      axis.title.x       = element_text(face = "bold", size = 13),
      axis.title.y       = element_text(face = "bold", size = 13),
      axis.text          = element_text(color = "black", size = 11),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", size = 0.3),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.position    = if (i == 1) c(0.6, 0.27) else "none",
      legend.background  = element_rect(fill = alpha("white", 0.5), color = NA),
      legend.justification = c(1, 1),
      legend.margin = margin(1, 1, 1, 1),
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.spacing.x = unit(0.2, "cm"),
      legend.spacing.y = unit(0.1, "cm")
    )
})

final_plot_pcr <- (plot_list[[1]] | plot_list[[2]] | plot_list[[3]]) /
  (plot_list[[4]] | plot_list[[5]] | plot_list[[6]])

ggsave("figures/FigS5.svg",
       plot = final_plot_pcr,
       width = 10, height = 6.5, dpi = 600)
