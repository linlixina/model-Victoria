#---------------------
# VIC epidemic context
#---------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(scales)
library(sf)
library(patchwork)
library(mgcv)

#Panel A
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



#Panel B
# Load the shapefile (replace 'path_to_shapefile' with the correct path to your downloaded file)
victoria_lga <- st_read("data/GDA94/vic_lga.shp")

# Updated list of LGAs to highlight
lga_list <- lga_list <- c(
  "Ararat Rural City", "Ballarat City", "Banyule City", "Bass Coast Shire", "Baw Baw Shire",
  "Bayside City", "Benalla Rural City", "Boroondara City", "Brimbank City", "Campaspe Shire",
  "Cardinia Shire", "Casey City", "Central Goldfields Shire", "Colac Otway Shire", "Corangamite Shire",
  "Darebin City", "East Gippsland Shire", "Frankston City", "Gannawarra Shire", "Glen Eira City",
  "Glenelg Shire", "Golden Plains Shire", "Greater Bendigo City", "Greater Dandenong City",
  "Greater Geelong City", "Greater Shepparton City", "Hepburn Shire", "Hobsons Bay City",
  "Hume City", "Indigo Shire", "Kingston City", "Knox City", "Latrobe City", "Macedon Ranges Shire",
  "Manningham City", "Maribyrnong City", "Maroondah City", "Melbourne City", "Melton City",
  "Mildura Rural City", "Mitchell Shire", "Moira Shire", "Monash City", "Moonee Valley City",
  "Moorabool Shire", "Merri-Bek City", "Mornington Peninsula Shire", "Mount Alexander Shire",
  "Moyne Shire", "Murrindindi Shire", "Nillumbik Shire", "Northern Grampians Shire",
  "Port Phillip City", "South Gippsland Shire", "Stonnington City", "Strathbogie Shire",
  "Surf Coast Shire", "Swan Hill Rural City", "Wangaratta Rural City", "Warrnambool City",
  "Wellington Shire", "Whitehorse City", "Whittlesea City", "Wodonga City", "Wyndham City",
  "Yarra City", "Yarra Ranges Shire"
)

victoria_lga$highlight <- ifelse(victoria_lga$LGA_NAME %in% lga_list, "Selected", "Other")

map_plot <- ggplot(data = victoria_lga) +
  geom_sf(aes(fill = highlight), color = "black") +
  scale_fill_manual(values = c("Selected" = "lightblue", "Other" = "gray90")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    legend.position  = c(0.85, 0.85),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 11),
    legend.key.size  = unit(0.5, "cm"),
    legend.background= element_rect(fill = "white", color = NA),
    axis.title       = element_blank(), 
    axis.text        = element_text(size = 11)
  )


#plotting!!!
combined <- (p + map_plot) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave("figures/Fig1.svg", combined, width = 11, height = 4.8, dpi = 600)
