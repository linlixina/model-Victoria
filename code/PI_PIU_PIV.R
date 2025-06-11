#-----------
# PI,PIU,PIV
#-----------
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(ggsci)

p <- 0.01
k <- 0.85
R0_list <- c(1.5, 4)
v_vals  <- seq(0, 1, length.out = 5000)

solve_Z <- function(v, R0, p, k, tol = 1e-6, max_iter = 1e4) {
  Z <- 0.5
  for (i in seq_len(max_iter)) {
    e1 <- exp(-R0 * (Z + p))
    e2 <- exp(- (1 - k) * R0 * (Z + p))
    F  <- Z - (1 - p) * (1 - (1 - v) * e1 - v * e2)
    dF <- 1 - (1 - p) * ((1 - v) * R0 * e1 + v * (1 - k) * R0 * e2)
    if (abs(dF) < tol) break
    Znew <- Z - F / dF
    if (abs(Znew - Z) < tol) break
    Z <- Znew
  }
  Z
}

find_roots <- function(x, ypp) {
  i <- which(ypp[-1] * ypp[-length(ypp)] < 0)[1]
  if (is.na(i)) return(NA_real_)
  x0 <- x[i]; x1 <- x[i+1]; y0 <- ypp[i]; y1 <- ypp[i+1]
  x0 - y0 * (x1 - x0) / (y1 - y0)
}


df <- expand.grid(v = v_vals, R0 = R0_list) %>%
  rowwise() %>%
  mutate(
    Z  = solve_Z(v, R0, p, k),
    e1 = exp(-R0*(Z + p)),
    e2 = exp(-(1 - k)*R0*(Z + p))
  ) %>%
  ungroup() %>%
  mutate(
    PI       = Z / (1 - p),
    PI_U     = 1 - e1,
    PI_V     = 1 - e2,
    Zp        = (1 - p)*(e1 - e2) / (1 - (1 - p)*R0*((1 - v)*e1 + v*(1 - k)*e2)),
    PI_p     = Zp / (1 - p),
    PI_U_p   = R0 * Zp * e1,
    PI_V_p   = (1 - k)*R0 * Zp * e2,
    Zpp       = {
      f_v <- (1 - p)*(e1 - e2)
      g_v <- 1 - (1 - p)*R0*((1 - v)*e1 + v*(1 - k)*e2)
      f_p <- (1 - p)*R0*Zp*(-e1 + (1 - k)*e2)
      g_p <- (1 - p)*R0*((1 + (1 - v)*R0*Zp)*e1 + (v*(1 - k)^2*R0*Zp - (1 - k))*e2)
      (f_p * g_v - f_v * g_p) / g_v^2
    },
    PI_pp     = Zpp / (1 - p),
    PI_U_pp   = R0 * e1 * (Zpp - R0 * Zp^2),
    PI_V_pp   = (1 - k)*R0 * e2 * (Zpp - (1 - k)*R0 * Zp^2)
  )

v_star <- df %>%
  group_by(R0) %>%
  summarise(
    v_star_PI   = find_roots(v, PI_pp),
    v_star_PI_U = find_roots(v, PI_U_pp),
    v_star_PI_V = find_roots(v, PI_V_pp)
  ) %>%
  ungroup() %>%
  pivot_longer(starts_with("v_star"), names_to = "nm", values_to = "vstar") %>%
  mutate(
    Function = recode(nm,
                      v_star_PI   = "bold(PI(v))",
                      v_star_PI_U = "bold(PI[U](v))",
                      v_star_PI_V = "bold(PI[V](v))"),
    Order = "2nd"
  )

df_long <- df %>%
  select(v, R0,
         PI,   PI_p,   PI_pp,
         PI_U, PI_U_p, PI_U_pp,
         PI_V, PI_V_p, PI_V_pp) %>%
  pivot_longer(-c(v, R0), names_to = "nm", values_to = "value") %>%
  mutate(
    Function = case_when(
      grepl("^PI_V", nm) ~ "bold(PI[V](v))",
      grepl("^PI_U", nm) ~ "bold(PI[U](v))",
      TRUE                ~ "bold(PI(v))"
    ),
    Order = case_when(
      grepl("_pp$", nm) ~ "2nd",
      grepl("_p$",  nm) ~ "1st",
      TRUE              ~ "0th"
    )
  )

ann <- expand.grid(
  Order    = c("0th","1st","2nd"),
  Function = c("bold(PI(v))","bold(PI[U](v))","bold(PI[V](v))")
) %>%
  arrange(Order, Function) %>%
  mutate(tag = LETTERS[1:9])

my_theme_final <- theme_minimal(base_size = 14) +
  theme(
    axis.title           = element_text(size = 13, face = "bold"),
    axis.text            = element_text(size = 11),
    
    strip.text.x         = element_text(size = 12, face = "bold", margin = margin(t = 4, b = 4)),
    strip.text.y         = element_text(size = 12, face = "bold", angle = 0, margin = margin(r = 4, l = 4)),
    strip.placement      = "outside",
    strip.background     = element_blank(),
    
    panel.spacing        = unit(1.2, "lines"),
    
    panel.background     = element_rect(fill = NA, colour = NA),
    panel.border         = element_rect(fill = NA, colour = "black", size = 0.6),
    
    panel.grid.minor     = element_blank(),
    panel.grid.major.x   = element_blank(),
    panel.grid.major.y   = element_line(color = "grey80", size = 0.4),
    
    axis.line            = element_line(color = "black", size = 0.6),
    axis.ticks           = element_line(color = "black", size = 0.4),
    
    legend.position      = c(0.22, 0.99),
    legend.justification = c(0,    1),
    legend.direction     = "vertical",
    legend.background    = element_blank(),
    legend.title         = element_text(size = 12, face = "bold"),
    legend.text          = element_text(size = 11),
    legend.key           = element_blank(),
    legend.key.size      = unit(0.8, "lines"),
    
    plot.margin          = margin(5, 5, 5, 5, "pt")
  )

p_plot <- ggplot(df_long,
                 aes(x = v, y = value,
                     colour = factor(R0),
                     group  = factor(R0)
                 )) +
  geom_line(size = 0.8) +
  geom_vline(data = v_star,
             aes(xintercept = vstar, colour = factor(R0)),
             linetype = "dotted", size = 0.6, alpha = 0.7,
             show.legend = FALSE) +
  geom_text(data = v_star,
            aes(x = vstar, y = 0, label = sprintf("%.3f", vstar)),
            vjust = -0.5, size = 3, show.legend = FALSE) +
  geom_text(data = ann,
            aes(x = -Inf, y = Inf, label = tag),
            hjust = -0.5, vjust = 1.5, size = 5, inherit.aes = FALSE) +
  facet_grid(
    Order ~ Function,
    scales = "free_y",
    switch = "y",
    labeller = labeller(Function = label_parsed)
  ) +
  scale_colour_npg(name = expression(R[0]), labels = c("1.5","4")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = "Overall vaccination level v", y = NULL) +
  my_theme_final +
  coord_cartesian(clip = "off")

print(p_plot)

ggsave("figures/FigS6.svg", p_plot,
       device = "svg", width = 8, height = 6, units = "in", dpi = 300)

