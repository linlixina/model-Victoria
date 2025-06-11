#-----------------------------------------------------
# A comparison of PI(0) estimators for R0 = 1.5 and 4.
#-----------------------------------------------------
#R0=1.5
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
library(rstan)
library(scales)
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
library(rstan)
library(scales)
library(patchwork)
library(mgcv)

R0_value <- 1.5
dp <- 0.01
k <- 0.85
v_vals <- seq(0, 1, length.out = 100)

solve_Z <- function(v, R0, p, k, tol = 1e-6, max_iter = 10000) {
  Z <- 0.5
  for (i in 1:max_iter) {
    exp1 <- exp(-R0 * (Z + p))
    exp2 <- exp(-(1 - k) * R0 * (Z + p))
    F <- Z - (1 - p) * (1 - (1 - v) * exp1 - v * exp2)
    dF <- 1 - (1 - p) * ((1 - v) * R0 * exp1 + v * (1 - k) * R0 * exp2)
    if (dF == 0) break
    Z_new <- Z - F / dF
    if (abs(Z_new - Z) < tol) break
    Z <- Z_new
  }
  Z
}

theory_df <- data.frame(v = v_vals) %>%
  mutate(
    Z = sapply(v, function(vv) solve_Z(vv, R0_value, dp, k)),
    SU_inf = (1 - dp) * (1 - v) * exp(-R0_value * (Z + dp)),
    SV_inf = (1 - dp) * v * exp(-(1 - k) * R0_value * (Z + dp)),
    infect_rate_U = ((1 - dp) * (1 - v) - SU_inf) / ((1 - dp) * (1 - v)),
    infect_rate_V = ((1 - dp) * v - SV_inf) / ((1 - dp) * v),
    PI_theory = infect_rate_U
  )

# Subset for Stan
subset_df <- theory_df %>% filter(v >= 0.4 & v <= 0.6)
wide_df <- subset_df %>% rename(PI_R0_1.5 = infect_rate_U)

# Conservative1: mean of the 20 data points
conservative1_value <- mean(head(subset_df$infect_rate_U, 20))

stan_code <- "
data {
  int<lower=1> N;
  vector[N] x;
  vector<lower=0, upper=1>[N] y;
}
parameters {
  real<lower=0, upper=1> L;
  real <upper=0> k;
  real x0;
  real<lower=0> phi;
}
transformed parameters {
  vector[N] mu;
  for (i in 1:N)
    mu[i] = L / (1 + exp(-k * (x[i] - x0)));
}
model {
  L ~ normal(0.15, 0.03);
  k ~ normal(0, 5);
  x0 ~ normal(0.5, 0.2);
  phi ~ gamma(0.01, 0.01);
  y ~ beta(mu * phi, (1 - mu) * phi);
}
"

df_stan <- list(
  N = nrow(wide_df),
  x = wide_df$v,
  y = wide_df$PI_R0_1.5
)
fit <- stan(model_code = stan_code,
            data = df_stan,
            iter = 10000,
            chains = 4,
            seed = 123,
            refresh = 0)
posterior <- rstan::extract(fit)
L_mean  <- mean(posterior$L)
k_mean  <- mean(posterior$k)
x0_mean <- mean(posterior$x0)

x_new <- seq(0, 1, length.out = 200)
draws <- length(posterior$L)
y_post <- sapply(1:draws, function(i) {
  L_i  <- posterior$L[i]
  k_i  <- posterior$k[i]
  x0_i <- posterior$x0[i]
  L_i / (1 + exp(-k_i * (x_new - x0_i)))
})
y_fit <- apply(y_post, 1, mean)
y_lower <- apply(y_post, 1, quantile, 0.025)
y_upper <- apply(y_post, 1, quantile, 0.975)

data_fit <- tibble(
  x = x_new,
  y_fit = y_fit,
  lower = y_lower,
  upper = y_upper
)
data_obs <- tibble(
  x = wide_df$v,
  y = wide_df$PI_R0_1.5
)
lin_mod <- lm(y ~ x, data = data_obs)
data_lin <- tibble(
  x    = x_new,
  y_lin = predict(lin_mod, newdata = data.frame(x = x_new))
)

#preparing_plot_R0=1.5
my_theme <- theme_classic(base_size = 14) +
  theme(
    panel.border      = element_rect(color = "black", fill = NA, size = 0.6),
    panel.grid        = element_blank(),
    axis.title        = element_text(size = 14, face = "bold"),
    axis.text         = element_text(size = 12),
    axis.line         = element_line(color = "black"),
    legend.position   = c(0.75, 0.83),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key        = element_blank(),
    legend.spacing.y  = unit(0.02, "cm"),
    legend.text       = element_text(size = 12),
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5)
  )

m= conservative1_value
n= data_fit$y_fit[80]
p_fit <- ggplot() +
  geom_point(data = data_obs,
             aes(x = x, y = y),
             shape = 21,
             color = "red",
             fill = "white",
             size = 3,
             stroke = 1.2,
             alpha = 1,
             show.legend = FALSE)+
  
  geom_line(data = theory_df,
            aes(x = v, y = PI_theory),
            color = "#4d4d4d", size = 1.0, linetype = "solid", show.legend = FALSE) +
  
  geom_line(data = data_fit,
            aes(x = x, y = y_fit),
            color = "#1b4f72", size = 1.0, show.legend = FALSE) +
  
  geom_ribbon(data = data_fit,
              aes(x = x, ymin = lower, ymax = upper),
              fill = "#1b4f72", alpha = 0.15, show.legend = FALSE) +
  
  geom_segment(aes(x = 0, xend = 0.4,
                   y = m, yend = m,
                   color = "CONSERVATIVE1", linetype = "CONSERVATIVE1"),
               size = 0.8) +
  
  geom_segment(aes(x = 0, xend = 0.4,
                   y = n, yend = n,
                   color = "CONSERVATIVE2", linetype = "CONSERVATIVE2"),
               size = 0.8) +
  
  geom_line(data = data_lin,
            aes(x = x, y = y_lin,
                color = "Linear fit", linetype = "Linear fit"),
            size = 0.9) +
  
  annotate("segment", x = 0.08, xend = 0.02, y = 0.25, yend = 0.19,
           arrow = arrow(length = unit(0.15, "inches")), color = "black") +
  annotate("text", x = 0.09, y = 0.25, label = "Logistic fit",
           hjust = 0, size = 4, color = "black") +
  
  annotate("segment", x = 0.18, xend = 0.12, y = 0.55, yend = 0.49,
           arrow = arrow(length = unit(0.15, "inches")), color = "black") +
  annotate("text", x = 0.19, y = 0.55, label = "SIR model",
           hjust = 0, size = 4, color = "black") +
  
  scale_color_manual(name = NULL,
                     values = c(
                       "Linear fit"     = "#d35400",
                       "CONSERVATIVE1"  = "black",
                       "CONSERVATIVE2"  = "black"
                     ),
                     breaks = c("CONSERVATIVE1", "CONSERVATIVE2", "Linear fit")) +
  
  scale_linetype_manual(name = NULL,
                        values = c(
                          "Linear fit"     = "dashed",
                          "CONSERVATIVE1"  = "dotted",
                          "CONSERVATIVE2"  = "dashed"
                        ),
                        breaks = c("CONSERVATIVE1", "CONSERVATIVE2", "Linear fit")) +
  
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 1)
  ) +
  
  labs(
    x = "Overall Vaccination Level v",
    y = expression(bold(PI[U](v))),
    title = expression(
      bold(R[0]) ~ bold("=") ~ bold(1.5)
    )
  ) +
  my_theme

p_A <- p_fit 






#R0=4.0
R0_value <- 4.0
dp <- 0.01
k <- 0.85
v_vals <- seq(0, 1, length.out = 100)

solve_Z <- function(v, R0, p, k, tol = 1e-6, max_iter = 10000) {
  Z <- 0.5
  for (i in 1:max_iter) {
    exp1 <- exp(-R0 * (Z + p))
    exp2 <- exp(-(1 - k) * R0 * (Z + p))
    F <- Z - (1 - p) * (1 - (1 - v) * exp1 - v * exp2)
    dF <- 1 - (1 - p) * ((1 - v) * R0 * exp1 + v * (1 - k) * R0 * exp2)
    if (dF == 0) break
    Z_new <- Z - F / dF
    if (abs(Z_new - Z) < tol) break
    Z <- Z_new
  }
  Z
}

theory_df <- data.frame(v = v_vals) %>%
  mutate(
    Z = sapply(v, function(vv) solve_Z(vv, R0_value, dp, k)),
    SU_inf = (1 - dp) * (1 - v) * exp(-R0_value * (Z + dp)),
    SV_inf = (1 - dp) * v * exp(-(1 - k) * R0_value * (Z + dp)),
    infect_rate_U = ((1 - dp) * (1 - v) - SU_inf) / ((1 - dp) * (1 - v)),
    infect_rate_V = ((1 - dp) * v - SV_inf) / ((1 - dp) * v),
    PI_theory = infect_rate_U
  )

subset_df <- theory_df %>% filter(v >= 0.4 & v <= 0.6)
wide_df <- subset_df %>% rename(PI_R0_1.5 = infect_rate_U)

conservative1_value <- mean(head(subset_df$infect_rate_U, 20))

stan_code <- "
data {
  int<lower=1> N;
  vector[N] x;
  vector<lower=0, upper=1>[N] y;
}
parameters {
  real<lower=0, upper=1> L;
  real<upper=0> k;
  real x0;
  real<lower=0> phi;
}
transformed parameters {
  vector[N] mu;
  for (i in 1:N)
    mu[i] = L / (1 + exp(-k * (x[i] - x0)));
}
model {
  L ~ normal(0.94, 0.03);
  k ~ normal(0, 5);
  x0 ~ normal(0.5, 0.2);
  phi ~ gamma(0.01, 0.01);
  y ~ beta(mu * phi, (1 - mu) * phi);
}
"

df_stan <- list(
  N = nrow(wide_df),
  x = wide_df$v,
  y = wide_df$PI_R0_1.5
)
fit <- stan(model_code = stan_code,
            data = df_stan,
            iter = 10000,
            chains = 4,
            seed = 123,
            refresh = 0)
posterior <- rstan::extract(fit)
L_mean  <- mean(posterior$L)
k_mean  <- mean(posterior$k)
x0_mean <- mean(posterior$x0)

x_new <- seq(0, 1, length.out = 200)
draws <- length(posterior$L)
y_post <- sapply(1:draws, function(i) {
  L_i  <- posterior$L[i]
  k_i  <- posterior$k[i]
  x0_i <- posterior$x0[i]
  L_i / (1 + exp(-k_i * (x_new - x0_i)))
})
y_fit <- apply(y_post, 1, mean)
y_lower <- apply(y_post, 1, quantile, 0.025)
y_upper <- apply(y_post, 1, quantile, 0.975)

data_fit <- tibble(
  x = x_new,
  y_fit = y_fit,
  lower = y_lower,
  upper = y_upper
)
data_obs <- tibble(
  x = wide_df$v,
  y = wide_df$PI_R0_1.5
)
lin_mod <- lm(y ~ x, data = data_obs)
data_lin <- tibble(
  x    = x_new,
  y_lin = predict(lin_mod, newdata = data.frame(x = x_new))
)


#preparing_plot_R0=4.0

my_theme <- theme_classic(base_size = 14) +
  theme(
    panel.border      = element_rect(color = "black", fill = NA, size = 0.6),
    panel.grid        = element_blank(),
    axis.title        = element_text(size = 14, face = "bold"),
    axis.text         = element_text(size = 12),
    axis.line         = element_line(color = "black"),
    legend.position   = c(0.25, 0.2),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key        = element_blank(),
    legend.spacing.y  = unit(0.02, "cm"),
    legend.text       = element_text(size = 12),
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5)
  )

p_fit <- ggplot() +
  geom_point(data = data_obs,
             aes(x = x, y = y),
             shape = 21,
             color = "red",
             fill = "white",
             size = 3,
             stroke = 1.2,
             alpha = 1,
             show.legend = FALSE)+
  
  geom_line(data = theory_df,
            aes(x = v, y = PI_theory),
            color = "#4d4d4d", size = 1.0, linetype = "solid", show.legend = FALSE) +
  
  geom_line(data = data_fit,
            aes(x = x, y = y_fit),
            color = "#1b4f72", size = 1.0, show.legend = FALSE) +
  
  geom_ribbon(data = data_fit,
              aes(x = x, ymin = lower, ymax = upper),
              fill = "#1b4f72", alpha = 0.15, show.legend = FALSE) +
  
  geom_segment(aes(x = 0, xend = 0.4,
                   y = conservative1_value, yend = conservative1_value,
                   color = "CONSERVATIVE1", linetype = "CONSERVATIVE1"),
               size = 0.8) +
  
  geom_segment(aes(x = 0, xend = 0.4,
                   y = data_fit$y_fit[80], yend = data_fit$y_fit[80],
                   color = "CONSERVATIVE2", linetype = "CONSERVATIVE2"),
               size = 0.8) +
  
  geom_line(data = data_lin,
            aes(x = x, y = y_lin,
                color = "Linear fit", linetype = "Linear fit"),
            size = 0.9) +
  
  annotate("segment", x = 0.1, xend = 0.02, y = 0.88, yend = 0.9579,
           arrow = arrow(length = unit(0.15, "inches")), color = "black") +
  annotate("text", x = 0.1, y = 0.87, label = "Logistic fit",
           hjust = 0, size = 4, color = "black") +
  
  annotate("segment", x = 0.17, xend = 0.08, y = 1.04, yend = 0.972,
           arrow = arrow(length = unit(0.15, "inches")), color = "black") +
  annotate("text", x = 0.17, y = 1.05, label = "SIR model",
           hjust = 0, size = 4, color = "black") +
  
  scale_color_manual(name = NULL,
                     values = c(
                       "Linear fit"     = "#d35400",
                       "CONSERVATIVE1"  = "black",
                       "CONSERVATIVE2"  = "black"
                     ),
                     breaks = c("CONSERVATIVE1", "CONSERVATIVE2", "Linear fit")) +
  
  scale_linetype_manual(name = NULL,
                        values = c(
                          "Linear fit"     = "dashed",
                          "CONSERVATIVE1"  = "dotted",
                          "CONSERVATIVE2"  = "dashed"
                        ),
                        breaks = c("CONSERVATIVE1", "CONSERVATIVE2", "Linear fit")) +
  
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 1)
  ) +
  
  labs(
    x = "Overall Vaccination Level v",
    y = expression(bold(PI[U](v))),
    title = expression(
      bold(R[0]) ~ bold("=") ~ bold(4.0)
    )
  ) +
  my_theme

p_B <- p_fit 





#plotting!!!
combined <- (p_A + p_B) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave("figures/Fig3.svg", combined, width = 11, height = 4.8, dpi = 600)

