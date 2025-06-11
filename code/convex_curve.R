#--------------------------------------------------------------
# Schematic illustration using a convex downward function PI(v)
#--------------------------------------------------------------
library(ggplot2)
library(grid)

xmin <- -0.1; xmax <- 0.75
ymin <- -0.0109; ymax <- 0.06

L <- 0.0592
k <- -13.39
x0 <- 0.42
f <- function(x) L / (1 + exp(-k * (x - x0)))

v1  <- 0.45; vm  <- 0.55; v2 <- 0.65
A1  <- c(v1,  f(v1));  Am <- c(vm, f(vm));  A2 <- c(v2, f(v2))

df_curve  <- data.frame(x = seq(0+0.22, xmax-0.04, length.out = 200))
df_curve$y <- f(df_curve$x)

df_vline  <- data.frame(x = c(v1, vm, v2), yend = c(A1[2], Am[2], A2[2]))
df_hline  <- data.frame(x = c(A1[1], Am[1], A2[1]), y = c(A1[2], Am[2], A2[2]))
df_points <- data.frame(x = c(A1[1], Am[1], A2[1]), y = c(A1[2], Am[2], A2[2]))

df_xlbl <- data.frame(
  x     = c(v1, vm, v2),
  y     = rep(-0.002, 3),
  label = c("bold(v[1])", "bold(v[m])", "bold(v[2])")
)

df_ylbl <- data.frame(
  x     = c(0.2-0.02,0.2-0.015,0.2-0.02),
  y     = c(A1[2], Am[2], A2[2]),
  label = c("bold(PI(v[1]))",
            "bold(PI(v[m]))",
            "bold(PI(v[2]))")
)

v_star <- 0.42
y_star <- f(v_star)
df_star <- data.frame(x = v_star, y = y_star)

df_vstar_vline <- data.frame(x = v_star, yend = y_star)
df_vstar_xlbl <- data.frame(
  x     = v_star,
  y     = -0.002+0.00052,
  label = "v*"
)

df_vstar_ylbl <- data.frame(
  x     = 0.2 - 0.015,
  y     = y_star,
  label = "bold(PI(v^\"*\"))"
)

p <- ggplot() +
  geom_blank(data = data.frame(x = c(xmin, xmax), y = c(ymin, ymax)),
             aes(x, y)) +
  
  geom_segment(aes(x = xmin+0.28, y = 0, xend = xmax, yend = 0),
               arrow = arrow(length = unit(0.07, "inches"), ends = "last"),
               size = 1.2) +
  
  geom_segment(aes(x = xmin+0.3, y = ymin+0.008, xend = xmin+0.3, yend = ymax),
               arrow = arrow(length = unit(0.07, "inches"), ends = "last"),
               size = 1.2) +
  
  geom_line(data = df_curve, aes(x, y), size = 1.2) +
  
  geom_segment(data = df_vline,
               aes(x = x, xend = x, y = 0, yend = yend),
               linetype = "dashed", color = "gray50") +
  
  geom_segment(data = df_hline,
               aes(x = xmin+0.3, xend = x, y = y, yend = y),
               linetype = "dashed", color = "gray50") +
  
  geom_segment(data = df_vstar_vline,
               aes(x = x, xend = x, y = 0, yend = yend),
               linetype = "dashed", color = "gray50") +
  
  
  geom_segment(aes(x = A1[1], y = A1[2], xend = Am[1], yend = Am[2]),
               arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
               color = "steelblue4", size = 1.8) +
  
  geom_segment(aes(x = A2[1], y = A2[2], xend = Am[1], yend = Am[2]),
               arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
               color = "firebrick", size = 1.8) +
  
  geom_point(data = df_points, aes(x, y), size = 4) +
  geom_point(data = df_star, aes(x, y), size = 4, color = "firebrick") +
  
  geom_text(data = df_xlbl, aes(x = x, y = y, label = label),
            parse = TRUE, vjust = 1, size = 5, fontface = "bold") +
  
  geom_text(data = df_ylbl, aes(x = x, y = y, label = label),
            parse = TRUE, hjust = 1, size = 5, fontface = "bold") +
  
  geom_text(data = df_vstar_xlbl, aes(x = x, y = y, label = label),
            vjust = 1, size = 5, fontface = "bold") +
  
  
  theme_void() +
  coord_cartesian(xlim = c(xmin, xmax),
                  ylim = c(ymin, ymax),
                  expand = FALSE,
                  clip = "off") +
  theme(plot.margin = unit(c(5, 5, -7.5, -45), "mm"))

ggsave("figures/Fig4.svg", p,
       width = 5.6, height = 4.4, units = "in", device = "svg")

