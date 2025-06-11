#---------------------------------------------------
# Sensitivity analysis of Poverty RATE and SEX RATIO
#---------------------------------------------------

#PI plots with Poverty RATE
library(rstan)
library(readxl)
library(ggplot2)
library(dplyr)
library(ggnewscale)
library(scales)
library(patchwork)

fit_lines_all <- data.frame()
points_all <- data.frame()

results_list <- readRDS("stanfits/all_results_list.rds")

df <- read_excel("data/Age_LGA.xlsx")
x <- df$overall_Vacc_Level_average
columns <- names(results_list)

fit_lines_all <- do.call(rbind, lapply(columns, function(col) {
  post  <- results_list[[col]]$posterior
  x_new <- seq(min(x), max(x), length.out = 200)

  y_mat <- sapply(seq_along(post$L), function(i) {
    L_i  <- post$L[i]
    k_i  <- post$k[i]
    x0_i <- post$x0[i]
    L_i / (1 + exp(-k_i * (x_new - x0_i)))
  })
  
  data.frame(
    x      = x_new,
    y_fit  = rowMeans(y_mat),
    lower  = apply(y_mat, 1, quantile, 0.025),
    upper  = apply(y_mat, 1, quantile, 0.975),
    group  = col,
    row.names = NULL
  )
}))

for (col in columns) {
  cat("Processing:", col, "\n")
  y <- df[[col]]
  y[y == 0] <- 1e-7 
  
  stan_data <- list(
    N = length(x),
    x = x,
    y = y
  )
  
  # Collect observed point data; include poverty_rate for color mapping
  points_all <- rbind(points_all, data.frame(
    x = x,
    y = y,
    group = col,
    poverty_rate = df$Poverty_Rate
  ))
}

ordered_levels <- columns

y_label_settings <- setNames(
  list(
    expression(bold(PI['<12y'])),
    expression(bold(PI['12-15y'])),
    expression(bold(PI['16-49y, V'])),
    expression(bold(PI['16-49y, U'])),
    expression(bold(PI['50+y, V'])),
    expression(bold(PI['50+y, U']))
  ),
  columns
)

y_max_table <- data.frame(group = ordered_levels,
                          y_max_fixed = c(0.05, 0.05, 0.012, 0.08, 0.012, 0.06))
group_labels <- setNames(LETTERS[1:6], ordered_levels)

fit_lines_all <- left_join(fit_lines_all, y_max_table, by = "group")
points_all    <- left_join(points_all,    y_max_table, by = "group")

fit_lines_all$group <- factor(fit_lines_all$group, levels = ordered_levels)
points_all$group    <- factor(points_all$group, levels = ordered_levels)


#plotting
plot_list <- lapply(seq_along(ordered_levels), function(i) {
  grp <- ordered_levels[i]
  df_fit <- fit_lines_all %>% filter(group == grp)
  df_pts <- points_all %>% filter(group == grp) %>% 
    mutate(poverty_rate = df$Poverty_Rate[match(x, df$overall_Vacc_Level_average)])
  y_max  <- unique(df_fit$y_max_fixed)
  
  p <- ggplot() +
    geom_blank(aes(x = x, y = y_max), data = df_fit) +
    geom_ribbon(data = df_fit, aes(x = x, ymin = lower, ymax = upper, fill = "PI (95% BCI)"),
                alpha = 0.3) +
    geom_line(data = df_fit, aes(x = x, y = y_fit, linetype = "PI (mean)"),
              color = "black", linewidth = 1) +
    geom_point(data = df_pts,
               aes(x = x, y = y, color = poverty_rate, shape = "Observed data"),
               size = 1.5, alpha = 0.8) +
    
    scale_color_gradientn(
      name = "Poverty Rate",
      colors = c("#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027"),
      limits = c(7, 21),
      breaks = seq(7, 21, length.out = 5),
      guide = if(i == 2) "colorbar" else "none" 
    ) +
    scale_fill_manual(NULL, values = c("PI (95% BCI)" = "#B0C4DE")) +
    scale_shape_manual(NULL, values = c("Observed data" = 16)) +
    scale_linetype_manual(NULL, values = c("PI (mean)" = 1)) +
    scale_x_continuous(limits = c(0.4, 0.65),
                       breaks = seq(0.4, 0.65, 0.05),
                       oob = oob_squish) +
    scale_y_continuous(breaks = pretty_breaks(n = 5)) +
    labs(
      title = group_labels[grp],
      x = "Overall Vaccination Level v",
      y = y_label_settings[[grp]]
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = -0.1),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.1, "cm"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.3, "cm")
    )
  
  if (i == 1) {
    p <- p +
      theme(
        legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.5), color = NA),
        legend.margin = margin(1, 1, 1, 1),
        legend.spacing.x = unit(0.2, "cm"),
        legend.spacing.y = unit(0.1, "cm")
      ) +
      guides(
        shape = guide_legend(order = 1, override.aes = list(
          color = "#7F3B08",
          size = 2
        )),
        linetype = guide_legend(order = 2),
        fill = guide_legend(order = 3)
      )
  } else if (i == 2) {
    p <- p +
      theme(
        legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.5), color = NA),
        legend.margin = margin(1, 1, 1, 1)
      ) +
      guides(
        color = guide_colorbar(title = "Poverty rate (%)"),
        fill = "none",
        linetype = "none",
        shape = "none"
      )
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
})

final_plot <- (plot_list[[1]] | plot_list[[2]] | plot_list[[3]]) /
  (plot_list[[4]] | plot_list[[5]] | plot_list[[6]])

print(final_plot)
ggsave("figures/FigS1.svg", plot = final_plot,
       width = 10, height = 6.5, dpi = 600)



#PI plots with SEX RATIO
fit_lines_all <- data.frame()
points_all <- data.frame()

results_list <- readRDS("stanfits/all_results_list.rds")
df <- read_excel("data/Age_LGA.xlsx")
x <- df$overall_Vacc_Level_average
columns <- names(results_list)

fit_lines_all <- do.call(rbind, lapply(columns, function(col) {
  post  <- results_list[[col]]$posterior
  x_new <- seq(min(x), max(x), length.out = 200)
  
  y_mat <- sapply(seq_along(post$L), function(i) {
    L_i  <- post$L[i]
    k_i  <- post$k[i]
    x0_i <- post$x0[i]
    L_i / (1 + exp(-k_i * (x_new - x0_i)))
  })
  
  data.frame(
    x      = x_new,
    y_fit  = rowMeans(y_mat),
    lower  = apply(y_mat, 1, quantile, 0.025),
    upper  = apply(y_mat, 1, quantile, 0.975),
    group  = col,
    row.names = NULL
  )
}))


for (col in columns) {
  cat("Processing:", col, "\n")
  y <- df[[col]]
  y[y == 0] <- 1e-7
  
  stan_data <- list(
    N = length(x),
    x = x,
    y = y
  )
  # Collect observed point data; include Sex_ratio for color mapping
  points_all <- rbind(points_all, data.frame(
    x = x,
    y = y,
    group = col,
    sex_ratio = df$`Sex ratio (males per 100 females)`
  ))
}

ordered_levels <- columns
y_label_settings <- setNames(
  list(
    expression(bold(PI['<12y'])),
    expression(bold(PI['12-15y'])),
    expression(bold(PI['16-49y, V'])),
    expression(bold(PI['16-49y, U'])),
    expression(bold(PI['50+y, V'])),
    expression(bold(PI['50+y, U']))
  ),
  columns
)

y_max_table <- data.frame(group = ordered_levels,
                          y_max_fixed = c(0.05, 0.05, 0.012, 0.08, 0.012, 0.06))
group_labels <- setNames(LETTERS[1:6], ordered_levels)

fit_lines_all <- left_join(fit_lines_all, y_max_table, by = "group")
points_all    <- left_join(points_all,    y_max_table, by = "group")

fit_lines_all$group <- factor(fit_lines_all$group, levels = ordered_levels)
points_all$group    <- factor(points_all$group, levels = ordered_levels)

plot_list <- lapply(seq_along(ordered_levels), function(i) {
  grp <- ordered_levels[i]
  df_fit <- fit_lines_all %>% filter(group == grp)
  df_pts <- points_all %>% filter(group == grp) %>% 
    mutate(sex_ratio = df$`Sex ratio (males per 100 females)`[match(x, df$overall_Vacc_Level_average)])
  y_max  <- unique(df_fit$y_max_fixed)
  
  p <- ggplot() +
    geom_blank(aes(x = x, y = y_max), data = df_fit) +
    geom_ribbon(data = df_fit, aes(x = x, ymin = lower, ymax = upper, fill = "PI (95% BCI)"),
                alpha = 0.3) +
    geom_line(data = df_fit, aes(x = x, y = y_fit, linetype = "PI (mean)"),
              color = "black", linewidth = 1) +
    geom_point(data = df_pts,
               aes(x = x, y = y, color = sex_ratio, shape = "Observed data"),
               size = 1.5, alpha = 0.8) +
    
    scale_color_gradientn(
      name = "Sex Ratio",
      colors = c("#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027"),
      limits = c(91, 106),
      breaks = seq(91, 106, length.out = 5),
      guide = if(i == 2) "colorbar" else "none"
    ) +
    scale_fill_manual(NULL, values = c("PI (95% BCI)" = "#B0C4DE")) +
    scale_shape_manual(NULL, values = c("Observed data" = 16)) +
    scale_linetype_manual(NULL, values = c("PI (mean)" = 1)) +
    scale_x_continuous(limits = c(0.4, 0.65),
                       breaks = seq(0.4, 0.65, 0.05),
                       oob = oob_squish) +
    scale_y_continuous(breaks = pretty_breaks(n = 5)) +
    labs(
      title = group_labels[grp],
      x = "Overall Vaccination Level v",
      y = y_label_settings[[grp]]
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = -0.1),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.1, "cm"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.3, "cm")
    )
  
  if (i == 1) {
    p <- p +
      theme(
        legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.5), color = NA),
        legend.margin = margin(1, 1, 1, 1),
        legend.spacing.x = unit(0.2, "cm"),
        legend.spacing.y = unit(0.1, "cm")
      ) +
      guides(
        shape = guide_legend(order = 1, override.aes = list(
          color = "#7F3B08",
          size = 2
        )),
        linetype = guide_legend(order = 2),
        fill = guide_legend(order = 3)
      )
  } else if (i == 2) {
    p <- p +
      theme(
        legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.5), color = NA),
        legend.margin = margin(1, 1, 1, 1)
      ) +
      guides(
        color = guide_colorbar(title = "Sex ratio (%)"),
        fill = "none",
        linetype = "none",
        shape = "none"
      )
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
})

final_plot <- (plot_list[[1]] | plot_list[[2]] | plot_list[[3]]) /
  (plot_list[[4]] | plot_list[[5]] | plot_list[[6]])

print(final_plot)
ggsave("figures/FigS2.svg", plot = final_plot,
       width = 10, height = 6.5, dpi = 600)

