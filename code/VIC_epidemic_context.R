#---------------------
# VIC epidemic context
#---------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(scales)

file_path   <- "data/VIC_epidemic.xlsx"
sheet1_name <- "Sheet1"
sheet2_name <- "Sheet2"

df1_raw <- read_excel(file_path, sheet = sheet1_name)
df2_raw <- read_excel(file_path, sheet = sheet2_name)

time1     <- names(df1_raw)[1]
cases_col <- names(df1_raw)[2]

time2     <- names(df2_raw)[1]
m1_col    <- names(df2_raw)[2]
m2_col    <- names(df2_raw)[3]
m3_col    <- names(df2_raw)[4]

#df1: Cases (2020-2021)
df1 <- df1_raw %>%
  mutate(
    date  = as.Date(.data[[time1]], format = "%d %b %y"),
    cases = as.numeric(gsub(",", "", as.character(.data[[cases_col]])))
  ) %>%
  filter(date >= as.Date("2020-01-01") &
           date <= as.Date("2021-12-31"))

#df2: Vaccination coverage
df2 <- df2_raw %>%
  mutate(
    date                       = as.Date(.data[[time2]], format = "%d %b %y"),
    `Second dose (overall)`    = as.numeric(gsub(",", "", as.character(.data[[m1_col]]))) * 100,
    `Second dose (16-49 years)`= as.numeric(gsub(",", "", as.character(.data[[m2_col]]))) * 100,
    `Second dose (50+ years)`  = as.numeric(gsub(",", "", as.character(.data[[m3_col]]))) * 100
  )

max_cases <- max(df1$cases, na.rm = TRUE)
max_vax   <- max(
  df2$`Second dose (overall)`,
  df2$`Second dose (16-49 years)`,
  df2$`Second dose (50+ years)`,
  na.rm = TRUE
)
ratio <- max_vax / max_cases

my_theme <- theme_classic(base_size = 14) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position  = c(0.3, 0.82),
    legend.spacing   = unit(0.02, "cm"),
    legend.margin    = margin(0, 0, 0, 0),
    legend.key.size  = unit(0.5, "cm"),
    legend.background= element_rect(fill = "white", color = NA),
    legend.key       = element_rect(fill = "white", color = NA),
    axis.title       = element_text(size = 12, face = "bold"),
    axis.text        = element_text(size = 11),
    plot.title       = element_text(size = 14, face = "bold", hjust = 0.5)
  )

p <- ggplot() +
  geom_rect(
    inherit.aes = FALSE,
    aes(
      xmin = as.Date("2021-09-26"),
      xmax = as.Date("2021-11-21"),
      ymin = -Inf, ymax = Inf
    ),
    fill = "grey80", alpha = 0.4
  ) +
  
  geom_line(
    data = df1,
    aes(x = date, y = cases, color = "Cases"),
    size = 0.8
  ) +

  geom_line(
    data = df2,
    aes(x = date, y = `Second dose (overall)` / ratio, color = "Second dose (overall)"),
    size = 0.8
  ) +
  
  geom_line(
    data = df2,
    aes(x = date, y = `Second dose (16-49 years)` / ratio, color = "Second dose (16-49 years)"),
    size = 0.8
  ) +
  
  geom_line(
    data = df2,
    aes(x = date, y = `Second dose (50+ years)` / ratio, color = "Second dose (50+ years)"),
    size = 0.8
  ) +

  scale_x_date(
    breaks = seq(as.Date("2020-01-01"), as.Date("2021-12-01"), by = "3 months"),
    labels = function(x) {
      ifelse(
        format(x, "%m") == "01",
        paste0(format(x, "%b"), "\n", format(x, "%Y")),
        format(x, "%b")
      )
    }
  ) +

  scale_y_continuous(
    name     = "Number of confirmed cases\ninfected with SARS-CoV-2",
    labels   = comma,
    sec.axis = sec_axis(
      ~ . * ratio,
      name   = "Vaccination coverage (%)",
      labels = comma
    )
  ) +

  scale_color_manual(
    "",
    values = c(
      "Cases"                     = "#4D4D4D",
      "Second dose (overall)"     = "#0072B2",
      "Second dose (16-49 years)" = "#009E73",
      "Second dose (50+ years)"   = "#D55E00"
    )
  ) +

  labs(x = NULL) +
  my_theme

ggsave(
  filename = "figures/Fig1.svg",
  plot     = p,
  width    = 6, height = 4,
  units    = "in",
  device   = "svg",
  dpi      = 600
)
