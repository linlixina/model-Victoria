#----------------------------------------
# Fitting with selected gamma1 and gamma2
#----------------------------------------

#fitting with selected sigma1 and sigma2
library(rstan)
library(readxl)
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)

df <- read_excel("data/Age_LGA.xlsx")
x <- df$overall_Vacc_Level_average

columns <- c(
  "Rate(Age<12)",  
  "Rate(12<=Age<16)",  
  "Rate(16=<V<50)", 
  "Rate(16=<UV<50)", 
  "Rate(V>=50)", 
  "Rate(UV>=50)"
)

stan_code <- "
data {
  int<lower=1> N;
  vector[N] x;
  vector<lower=0, upper=1>[N] y;
}
parameters {
  real<lower=0, upper=1> L;
  real k;
  real x0;
  real<lower=0> phi;
}
transformed parameters {
  vector[N] mu;
  for (i in 1:N)
    mu[i] = L / (1 + exp(-k * (x[i] - x0)));
}
model {
  L ~ normal(0.01, 0.003);
  k ~ normal(0, 5);
  x0 ~ normal(0.5, 0.2);
  phi ~ gamma(0.01, 0.01);
  y ~ beta(mu * phi, (1 - mu) * phi);
}
"

stan_code_1 <- "
data {
  int<lower=1> N;
  vector[N] x;
  vector<lower=0, upper=1>[N] y;
}
parameters {
  real<lower=0, upper=1> L;
  real k;
  real x0;
  real<lower=0> phi;
}
transformed parameters {
  vector[N] mu;
  for (i in 1:N)
    mu[i] = L / (1 + exp(-k * (x[i] - x0)));
}
model {
  L ~ normal(0.05, 0.015);
  k ~ normal(0, 5);
  x0 ~ normal(0.5, 0.2);
  phi ~ gamma(0.01, 0.01);
  y ~ beta(mu * phi, (1 - mu) * phi);
}
"

output_dir <- "stanfits"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

results_list <- list()

for (col in columns) {
  cat("Processing:", col, "...\n")
  
  y <- df[[col]]
  y[y == 0] <- 1e-7
  
  stan_data <- list(
    N = length(x),
    x = x,
    y = y
  )
  
  if (col %in% c("Rate(16=<V<50)", "Rate(V>=50)")) {
    fit <- stan(model_code = stan_code,
                data       = stan_data,
                iter       = 10000,
                chains     = 4,
                seed       = 123,
                refresh    = 0)
  } else {
    fit <- stan(model_code = stan_code_1,
                data       = stan_data,
                iter       = 10000,
                chains     = 4,
                seed       = 123,
                refresh    = 0)
  }
  
  posterior <- rstan::extract(fit)
  L_mean <- mean(posterior$L)
  k_mean <- mean(posterior$k)
  x0_mean <- mean(posterior$x0)
  x0_q <- quantile(posterior$x0, probs = c(0.025, 0.975))
  
  results_list[[col]] <- list(
    group = col,
    x = x,
    y = y,
    fit = fit,
    posterior = posterior,
    L_mean = L_mean,
    k_mean = k_mean,
    x0_mean = x0_mean,
    x0_q = x0_q
  )
}

saveRDS(results_list, file = file.path(output_dir, "all_results_list.rds"))



#-------------------------------------------------
# plotting PI fitted curve using posterior results
#-------------------------------------------------
#Read Stan model output (rds file)
results_list <- readRDS("stanfits/all_results_list.rds")

fit_lines_obs  <- data.frame()
points_all_obs <- data.frame()

for (name in names(results_list)) {
  res   <- results_list[[name]]
  x_new <- seq(min(x), max(x), length.out = 200)
  draws <- length(res$posterior$L)

  y_post <- sapply(1:draws, function(i) {
    L_i  <- res$posterior$L[i]
    k_i  <- res$posterior$k[i]
    x0_i <- res$posterior$x0[i]
    L_i / (1 + exp(-k_i * (x_new - x0_i)))
  })
  
  y_fit   <- apply(y_post, 1, mean)
  y_lower <- apply(y_post, 1, quantile, 0.025)
  y_upper <- apply(y_post, 1, quantile, 0.975)
  
  fit_lines_obs <- rbind(fit_lines_obs, data.frame(
    x      = x_new,
    y_fit  = y_fit,
    lower  = y_lower,
    upper  = y_upper,
    group  = name
  ))
  
  points_all_obs <- rbind(points_all_obs, data.frame(
    x      = res$x,
    y      = res$y,
    group  = name
  ))
}

fit_lines_full   <- data.frame()
points_all_full  <- data.frame()
x0_summary_all   <- data.frame()
x_new_full       <- seq(0, 1, length.out = 1000)

for (name in names(results_list)) {
  res   <- results_list[[name]]
  draws <- length(res$posterior$L)
  
  y_post <- sapply(1:draws, function(i) {
    L_i  <- res$posterior$L[i]
    k_i  <- res$posterior$k[i]
    x0_i <- res$posterior$x0[i]
    L_i / (1 + exp(-k_i * (x_new_full - x0_i)))
  })
  
  y_fit   <- apply(y_post, 1, mean)
  y_lower <- apply(y_post, 1, quantile, 0.025)
  y_upper <- apply(y_post, 1, quantile, 0.975)
  
  fit_lines_full <- rbind(fit_lines_full, data.frame(
    x      = x_new_full,
    y_fit  = y_fit,
    lower  = y_lower,
    upper  = y_upper,
    group  = name
  ))
  
  points_all_full <- rbind(points_all_full, data.frame(
    x      = res$x,
    y      = res$y,
    group  = name
  ))
  
  x0_summary_all <- rbind(x0_summary_all, data.frame(
    group     = name,
    x0_mean   = res$x0_mean,
    x0_lower  = res$x0_q[1],
    x0_upper  = res$x0_q[2]
  ))
}

ordered_levels <- columns

y_label_settings_pcr <- list(
  "Rate(Age<12)"     = expression(bold(PI['<12y'])),
  "Rate(12<=Age<16)" = expression(bold(PI['12-15y'])),
  "Rate(16=<V<50)"   = expression(bold(PI['16-49y, V'])),
  "Rate(16=<UV<50)"  = expression(bold(PI['16-49y, U'])),
  "Rate(V>=50)"      = expression(bold(PI['50+y, V'])),
  "Rate(UV>=50)"     = expression(bold(PI['50+y, U']))
)

y_max_table <- data.frame(
  group       = ordered_levels,
  y_max_fixed = c(0.05, 0.05, 0.012, 0.08, 0.012, 0.06)
)

group_labels <- setNames(LETTERS[1:6], ordered_levels)

fit_lines_obs  <- left_join(fit_lines_obs,  y_max_table, by = "group")
points_all_obs <- left_join(points_all_obs, y_max_table, by = "group")

fit_lines_obs$group  <- factor(fit_lines_obs$group,  levels = ordered_levels)
points_all_obs$group <- factor(points_all_obs$group, levels = ordered_levels)

plot_list <- lapply(seq_along(ordered_levels), function(i) {
  grp    <- ordered_levels[i]
  df_fit <- filter(fit_lines_obs,  group == grp)
  df_pts <- filter(points_all_obs, group == grp)
  y_max  <- unique(df_fit$y_max_fixed)
  
  p <- ggplot() +
    geom_blank(aes(x = x, y = y_max), data = df_fit) +
    geom_ribbon(
      data = df_fit,
      aes(x = x, ymin = lower, ymax = upper, fill = "PI (95% BCI)"),
      alpha = 0.3
    ) +
    geom_line(
      data = df_fit,
      aes(x = x, y = y_fit, color = "PI (mean)"),
      linewidth = 1
    ) +
    geom_point(
      data = df_pts,
      aes(x = x, y = y, shape = "Observed data"),
      color = "#7F3B08", size = 1.5, alpha = 0.8
    ) +
    scale_color_manual(NULL, values = c("PI (mean)" = "black")) +
    scale_fill_manual(NULL, values = c("PI (95% BCI)" = "#B0C4DE")) +
    scale_shape_manual(NULL, values = c("Observed data" = 16)) +
    scale_x_continuous(
      limits = c(0.4, 0.65),
      breaks = seq(0.4, 0.65, 0.05),
      oob = oob_squish
    ) +
    scale_y_continuous(breaks = pretty_breaks(n = 5)) +
    labs(
      title = group_labels[grp],
      x     = "Overall Vaccination Level v",
      y     = y_label_settings_pcr[[grp]]
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title         = element_text(size = 13, face = "bold", hjust = -0.1),
      axis.title.x       = element_text(size = 12, face = "bold"),
      axis.title.y       = element_text(size = 12, face = "bold"),
      axis.text.x        = element_text(size = 10),
      axis.text.y        = element_text(size = 10, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor    = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
      axis.line         = element_line(color = "black", linewidth = 0.5),
      axis.ticks        = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.1, "cm"),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.text       = element_text(size = 10),
      legend.key.size   = unit(0.3, "cm")
    )
  
  if (i == 1) {
    p <- p +
      theme(
        legend.position      = c(0.95, 0.95),
        legend.justification = c(1, 1),
        legend.background    = element_rect(fill = alpha("white", 0.5), color = NA),
        legend.margin        = margin(1, 1, 1, 1),
        legend.spacing.x     = unit(0.2, "cm"),
        legend.spacing.y     = unit(0.1, "cm")
      ) +
      guides(
        shape = guide_legend(order = 1,
                             override.aes = list(shape = 16, color = "#7F3B08", size = 2, stroke = 0)
        ),
        color = guide_legend(order = 2,
                             override.aes = list(linetype = 1, shape = NA, color = "black", size = 1)
        ),
        fill = guide_legend(order = 3,
                            override.aes = list(fill = "#B0C4DE", shape = NA, color = NA)
        )
      )
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
})

final_plot_pcr <- (plot_list[[1]] | plot_list[[2]] | plot_list[[3]]) /
  (plot_list[[4]] | plot_list[[5]] | plot_list[[6]])

ggsave(
  filename = "figures/Fig2.svg",
  plot     = final_plot_pcr,
  width    = 10,
  height   = 6.5,
  dpi      = 600
)




#-----------------------------------------------------------------
# plotting PI (v from 0 to 1) fitted curve using posterior results
#-----------------------------------------------------------------
ordered_levels_full <- ordered_levels
y_axis_settings <- list(
  "Rate(Age<12)"     = list(breaks = c(0, 0.025, 0.05, 0.075, 0.1), limits = c(0, 0.1)),
  "Rate(12<=Age<16)" = list(breaks = c(0, 0.025, 0.05, 0.075, 0.1), limits = c(0, 0.1)),
  "Rate(16=<V<50)"   = list(breaks = c(0, 0.005, 0.01, 0.015, 0.02), limits = c(0, 0.02)),
  "Rate(16=<UV<50)"  = list(breaks = c(0, 0.025, 0.05, 0.075, 0.1), limits = c(0, 0.1)),
  "Rate(V>=50)"      = list(breaks = c(0, 0.005, 0.01, 0.015, 0.02), limits = c(0, 0.02)),
  "Rate(UV>=50)"     = list(breaks = c(0, 0.025, 0.05, 0.075, 0.1), limits = c(0, 0.1))
)

y_label_settings_full <- list(
  "Rate(Age<12)"     = expression(bold(PI['<12y'])),
  "Rate(12<=Age<16)" = expression(bold(PI['12-15y'])),
  "Rate(16=<V<50)"   = expression(bold(PI['16-49y, V'])),
  "Rate(16=<UV<50)"  = expression(bold(PI['16-49y, U'])),
  "Rate(V>=50)"      = expression(bold(PI['50+y, V'])),
  "Rate(UV>=50)"     = expression(bold(PI['50+y, U']))
)

x0_summary_all$group   <- factor(x0_summary_all$group,   levels = ordered_levels_full)
fit_lines_full$group   <- factor(fit_lines_full$group,   levels = ordered_levels_full)
points_all_full$group  <- factor(points_all_full$group,  levels = ordered_levels_full)

plot_list_x0 <- lapply(seq_along(ordered_levels_full), function(i) {
  grp      <- ordered_levels_full[i]
  df_line  <- filter(fit_lines_full,  group == grp)
  df_point <- filter(points_all_full, group == grp)
  x0_info  <- filter(x0_summary_all,  group == grp)
  yset     <- y_axis_settings[[grp]]
  
  show_legend <- (i == 2)
  
  p <- ggplot() +
    geom_vline(
      data = x0_info,
      aes(xintercept = x0_mean, linetype = "x0 (mean)"),
      color   = "orange3", linewidth = 1, show.legend = FALSE
    ) +
    geom_hline(
      data = x0_info,
      aes(yintercept = x0_mean, linetype = "x0 (mean)"),
      color   = "orange3", linewidth = 1
    ) +
    geom_rect(
      data = x0_info,
      aes(xmin = x0_lower, xmax = x0_upper, alpha = "x0 (95% BCI)"),
      ymin   = -Inf, ymax = Inf, fill = "#FFD700"
    ) +
    geom_ribbon(
      data = df_line,
      aes(x = x, ymin = lower, ymax = upper, fill = "PI (95% BCI)"),
      alpha = 0.3
    ) +
    geom_line(
      data = df_line,
      aes(x = x, y = y_fit, color = "PI (mean)"),
      linewidth = 1.1
    ) +
    geom_point(
      data = df_point,
      aes(x = x, y = y, shape = "Observed data"),
      color = "#7F3B08", size = 1.5, alpha = 0.8
    ) +
    scale_color_manual(
      values = c("PI (mean)" = "black"),
      labels = c("PI (mean)")
    ) +
    scale_shape_manual(
      values = c("Observed data" = 16),
      labels = c("Observed data")
    ) +
    scale_linetype_manual(
      values = c("x0 (mean)" = "solid"),
      labels = c(expression(x[0] * " (mean)"))
    ) +
    scale_fill_manual(
      values = c("PI (95% BCI)" = "#B0C4DE"),
      labels = c("PI (95% BCI)")
    ) +
    scale_alpha_manual(
      values = c("x0 (95% BCI)" = 0.25),
      labels = c(expression(x[0] * " (95% BCI)"))
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(limits = yset$limits, breaks = yset$breaks) +
    labs(
      title = LETTERS[i],
      x     = "Overall Vaccination Level v",
      y     = y_label_settings_full[[grp]] 
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title         = element_text(size = 14, face = "bold", hjust = -0.1),
      axis.title.x       = element_text(size = 13, face = "bold"),
      axis.title.y       = element_text(size = 13, face = "bold"), 
      axis.text.x        = element_text(size = 11),
      axis.text.y        = element_text(size = 11, face = "bold"), 
      panel.grid.major.x = element_blank(),
      panel.grid.minor    = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", size = 0.3),
      axis.line         = element_line(color = "black", linewidth = 0.6),
      axis.ticks.x      = element_line(color = "black", size = 0.4),
      axis.ticks.y      = element_line(color = "black", size = 0.4),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.position   = if (show_legend) c(0.995, 0.995) else "none",
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = alpha("white", 0.5), color = NA),
      legend.margin     = margin(1, 1, 1, 1),
      legend.key        = element_blank(),
      legend.title      = element_blank(),
      legend.spacing.x  = unit(0.2, "cm"),
      legend.spacing.y  = unit(0.1, "cm")
    )
  
  if (show_legend) {
    p <- p +
      guides(
        shape    = guide_legend(order = 1),
        color    = guide_legend(order = 2),
        linetype = guide_legend(order = 3),
        fill     = guide_legend(order = 4),
        alpha    = guide_legend(order = 5)
      )
  }
  
  return(p)
})

final_plot_x0 <- (plot_list_x0[[1]] | plot_list_x0[[2]] | plot_list_x0[[3]]) /
  (plot_list_x0[[4]] | plot_list_x0[[5]] | plot_list_x0[[6]])

ggsave(
  filename = "figures/FigS3.svg",
  plot     = final_plot_x0,
  width    = 10,
  height   = 6.5,
  dpi      = 600
)

