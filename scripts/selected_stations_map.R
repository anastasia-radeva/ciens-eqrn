library(dplyr)
library(ggplot2)
library(maps)
source(file = paste0(getwd(), "/init_file.R"))


# 1) Get station metadata 
stations_all <- loc_data %>%
  transmute(
    station_id = as.character(station_id),
    lon = as.numeric(longitude),   # or "longitude"
    lat = as.numeric(latitude)    # or "latitude"
  ) %>%
  filter(is.finite(lon), is.finite(lat)) %>%
  distinct(station_id, .keep_all = TRUE)

sel_ids <- c("10020","10453","10731","10963")

stations_all <- stations_all %>%
  mutate(selected = station_id %in% sel_ids)

# 2) Germany outline
de_map <- map_data("world", region = "Germany")

# 3) Plot
ggplot() +
  geom_polygon(
    data = de_map,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey40",
    linewidth = 0.3
  ) +
  geom_point(
    data = stations_all %>% filter(!selected),
    aes(x = lon, y = lat),
    size = 1.2,
    alpha = 0.6
  ) +
  geom_point(
    data = stations_all %>% filter(selected),
    aes(x = lon, y = lat),
    color = "maroon",
    size = 2.8
  ) +
  coord_quickmap(xlim = c(5, 16), ylim = c(47, 56)) +
  labs(
    title = "CIENS Stations",
    subtitle = "Selected Stations Highlighted",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal(base_family = "Franklin") + 
  theme(
  plot.title = element_text(size = 12),
  plot.subtitle = element_text(size = 10),
  axis.title = element_text(size = 10),
  axis.text  = element_text(size = 10),
  )
