#-------------------------------------------------------
# Sensitivity analysis of prior sigma for the two groups
#-------------------------------------------------------

library(rstan)
library(readxl)
library(dplyr)
library(tibble)

stan_code <- "
data {
  int<lower=1> N;
  vector[N] x;
  vector<lower=0,upper=1>[N] y;
  real<lower=0> prior_mu_L;
  real<lower=0> prior_sigma_L;
}
parameters {
  real<lower=0,upper=1> L;
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
  L   ~ normal(prior_mu_L, prior_sigma_L);
  k   ~ normal(0, 5);
  x0  ~ normal(0.5, 0.2);
  phi ~ gamma(0.01, 0.01);
  y   ~ beta(mu * phi, (1 - mu) * phi);
}
generated quantities {
  real PCR0;
  PCR0 = L / (1 + exp(-k * (0 - x0)));
}
"

df <- read_excel("data/Age_LGA.xlsx")
x  <- df$overall_Vacc_Level_average

groups_code1 <- c("Rate(16=<V<50)", "Rate(V>=50)")      # prior_mu_L = 0.01
groups_code2 <- c("Rate(Age<12)", "Rate(12<=Age<16)",
                  "Rate(16=<UV<50)", "Rate(UV>=50)")   # prior_mu_L = 0.05

sigma1_list <- c(0.01, 0.006, 0.003, 0.002, 0.001, 0.0005)    
sigma2_list <- c(0.05, 0.03, 0.015, 0.01, 0.005, 0.0025)      

sens_results <- tibble()

for (sigma in sigma1_list) {
  
  if (sigma <= 0.001) {
    iters  <- 10000
    warmup <- 5000
  } else {
    iters  <- 4000
    warmup <- 1000
  }
  
  for (grp in groups_code1) {
    
    y <- df[[grp]]
    y[y == 0] <- 1e-7
    data_list <- list(
      N             = length(x),
      x             = x,
      y             = y,
      prior_mu_L    = 0.01,
      prior_sigma_L = sigma
    )
    
    fit <- stan(
      model_code = stan_code,
      data       = data_list,
      iter       = iters,
      warmup     = warmup,
      chains     = 4,
      seed       = 123,
      refresh    = 0
    )
    
    sumdf <- summary(fit, pars = c("L", "PCR0"))$summary %>%
      as.data.frame() %>%
      rownames_to_column("parameter") %>%
      filter(parameter %in% c("L", "PCR0")) %>%
      mutate(group          = grp,
             prior_mu_L     = 0.01,
             prior_sigma_L  = sigma)
    
    sens_results <- bind_rows(sens_results, sumdf)
  }
}

for (sigma in sigma2_list) {
  
  if (sigma <= 0.005) {
    iters  <- 10000
    warmup <- 5000
  } else {
    iters  <- 4000
    warmup <- 1000
  }
  
  for (grp in groups_code2) {
    
    y <- df[[grp]]
    y[y == 0] <- 1e-7
    data_list <- list(
      N             = length(x),
      x             = x,
      y             = y,
      prior_mu_L    = 0.05,
      prior_sigma_L = sigma
    )
    
    fit <- stan(
      model_code = stan_code,
      data       = data_list,
      iter       = iters,
      warmup     = warmup,
      chains     = 4,
      seed       = 123,
      refresh    = 0
    )
    
    sumdf <- summary(fit, pars = c("L", "PCR0"))$summary %>%
      as.data.frame() %>%
      rownames_to_column("parameter") %>%
      filter(parameter %in% c("L", "PCR0")) %>%
      mutate(group         = grp,
             prior_mu_L    = 0.05,
             prior_sigma_L = sigma)
    
    sens_results <- bind_rows(sens_results, sumdf)
  }
}

sens_results <- sens_results %>%
  rename(
    mean_post = `mean`,
    ci_l = `2.5%`,
    ci_u = `97.5%`
    
  )

write.csv(sens_results,
          file     = "results/sensitivity_results_sigmas.csv",
          row.names = FALSE)





# plotting
library(dplyr)
library(ggplot2)
library(ggrepel)

group_levels <- c(
  "Rate(Age<12)", "Rate(12<=Age<16)",
  "Rate(16=<V<50)", "Rate(16=<UV<50)",
  "Rate(V>=50)", "Rate(UV>=50)"
)

panel_labels <- c(
  "Rate(Age<12)"     = "A",
  "Rate(12<=Age<16)" = "B",
  "Rate(16=<V<50)"   = "C",
  "Rate(16=<UV<50)"  = "D",
  "Rate(V>=50)"      = "E",
  "Rate(UV>=50)"     = "F"
)

group_labels <- c(
  "Rate(Age<12)"     = "<12 years",
  "Rate(12<=Age<16)" = "12–15 years",
  "Rate(16=<V<50)"   = "16–49 years, FV",
  "Rate(16=<UV<50)"  = "16–49 years, UV",
  "Rate(V>=50)"      = "50+ years, FV",
  "Rate(UV>=50)"     = "50+ years, UV"
)

axis_limits <- data.frame(
  group = group_levels,
  x_min = c(0, 0, 0, 0, 0, 0),
  x_max = c(0.05, 0.05, 0.01, 0.05, 0.01, 0.05),
  y_max = c(0.2, 0.2, 0.04, 0.2, 0.04, 0.2)
)

# Load sensitivity results
summary_tbl <- read.csv("results/sensitivity_results_sigmas.csv",
                        stringsAsFactors = FALSE)

summary_tbl <- summary_tbl %>%
  mutate(
    ratio = (ci_u - mean_post) / mean_post,
    group = factor(group, levels = group_levels),
    label = panel_labels[group],
    y_min = 0
  ) %>%
  left_join(axis_limits, by = "group") %>%
  mutate(
    err_width = 0.02 * (x_max - x_min)
  )

summary_tbl$group <- factor(summary_tbl$group, levels = group_levels)

label_df <- summary_tbl %>%
  filter(parameter == "L") %>%
  group_by(group, label) %>%
  summarise(x = min(x_min), y = max(y_max), .groups = "drop") %>%
  mutate(group = factor(group, levels = group_levels))

p_L_annot <- summary_tbl %>%
  filter(parameter == "L") %>%
  ggplot(aes(x = prior_sigma_L, y = mean_post, color = factor(prior_mu_L))) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u, width = err_width), linewidth = 0.4) +
  geom_blank(aes(x = x_min)) +
  geom_blank(aes(x = x_max)) +
  geom_blank(aes(y = y_max)) +
  geom_blank(aes(y = y_min)) +
  geom_text(data = label_df,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            hjust = -0.3, vjust = 1.3,
            fontface = "bold", size = 5) +
  facet_wrap(~ group, scales = "free", ncol = 3,
             labeller = labeller(group = group_labels)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  labs(
    x = expression(bold("Prior SD of") ~ bold(sigma[L])),
    y = "Posterior estimate of L (mean ± 95% BCI)",
    color = expression("Prior mean" ~ mu[L])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(face = "bold"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

print(p_L_annot)

ggsave(
  filename = "figures/FigS4.svg",
  plot = p_L_annot,
  width = 10,
  height = 6,
  dpi = 600,
  device = "svg"
)
