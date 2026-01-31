## CIENS: Data Exploration
## 12 UTC forecasts

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(purrr)
library(forcats)

# ------------------------------
#### Initialization ####
# ------------------------------

# Path to Github repository functions
git_path <- "~/"

# Set working directory
setwd(git_path)

# Source R file for initializing
source(file = paste0(getwd(), "/init_file.R"))


#### Get data: Construct data.frame from netcdf files according to your needs ####

# Define period of observation/forecast dates that are supposed to be in the created data.frame.
# In this example: forecasts/observation dates in the 2023.
tm_vec <- tm_vec[year(tm_vec) >= 2010 & year(tm_vec) <= 2023]

# Define observational and meteorological variables of interest, 
# here wind gust forecasts and observations.
obs_vars <- c("wind_speed_of_gust", "wind_speed")

# Generate data.frame df for all initializations
file_name <- file.path("~/forecasts/all_forecast_wind_gust")
if (file.exists(file_name)) {
  load(file = "~/forecasts/all_forecast_wind_gust")
  
} else {
  df <- bind_rows(lapply(tm_vec, function(x) get_init(tm = x,
                                                      dir_path = data_path,
                                                      obs_vars = obs_vars,
                                                      step_vec = c(0:21),  
                                                      ens_vec = ens_vec)))
  save(df ,file = "~/all_forecast_wind_gust")
}

# ------------------------------
#### Functions
# ------------------------------

#### Function which computes quantiles ####
station_quantile_matrix <- function(df, var = "wind_speed_of_gust",
                                    by = "location",
                                    probs = seq(0.01, 0.99, by = 0.01)) {
  stopifnot(var %in% names(df), by %in% names(df))
  x <- df[!is.na(df[[var]]), c(by, var)]
  # compute quantiles per station
  q_list <- tapply(x[[var]], x[[by]], function(v) as.numeric(quantile(v, probs = probs, na.rm = TRUE)))
  # drop stations with all-NA or too few data
  keep <- vapply(q_list, function(q) all(is.finite(q)), logical(1))
  q_list <- q_list[keep]
  Q <- do.call(rbind, q_list)           # stations x quantiles
  attr(Q, "stations") <- names(q_list)  # keep names
  Q
}

#### Function which chooses 4 stations that maximize the sum of pairwise distances ####
# Q:   station x quantile matrix
# exact: use exact combinatorial search or greedy heuristic
# fixed_stations: optional character vector, up to 2 station names to force-include
# weights: optional numeric vector of length ncol(Q) with non-negative weights
#          (e.g. larger weights for extreme quantiles)
top4_most_heterogeneous <- function(Q,
                                    exact = TRUE,
                                    fixed_stations = NULL,
                                    weights = NULL) {
  st <- attr(Q, "stations")
  if (is.null(st)) st <- rownames(Q)
  if (is.null(st)) stop("Q must have station names via attr() or rownames().")
  
  # --- handle weights ------------------------------------------------------
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      if (nrow(weights) == 1L || ncol(weights) == ncol(Q)) {
        weights <- as.numeric(weights)
      } else {
        stop("If 'weights' is a matrix, its length must match ncol(Q).")
      }
    }
    if (length(weights) != ncol(Q))
      stop("Length of 'weights' must equal ncol(Q).")
    if (any(weights < 0))
      stop("'weights' must be non-negative.")
    if (all(weights == 0))
      stop("All weights are zero; at least one weight must be > 0.")
    
    Wsqrt <- sqrt(weights)
    Q_use <- sweep(Q, 2, Wsqrt, `*`)
  } else {
    Q_use <- Q
  }
  
  # distance matrix
  D <- as.matrix(dist(Q_use, method = "euclidean"))
  diag(D) <- 0
  
  # --- handle fixed stations ------------------------------------------------
  if (!is.null(fixed_stations)) {
    fixed_stations <- intersect(fixed_stations, st)
    if (length(fixed_stations) > 3)
      stop("You can specify at most 3 fixed stations when selecting 4.")
    if (length(fixed_stations) == 0)
      fixed_stations <- NULL
  }
  
  fixed_idx <- which(st %in% fixed_stations)
  
  # --- exact search ---------------------------------------------------------
  if (exact) {
    cmb <- combn(seq_along(st), 4)
    
    if (length(fixed_idx) > 0) {
      keep <- apply(cmb, 2, function(x) all(fixed_idx %in% x))
      cmb <- cmb[, keep, drop = FALSE]
      if (ncol(cmb) == 0)
        stop("No quadruplets contain all specified 'fixed_stations'.")
    }
    
    # compute all six pairwise distances in the quadruple
    pair_idx <- combn(1:4, 2)  # all 6 pairs of (1,2,3,4)
    score <- apply(cmb, 2, function(g) {
      sum(D[g[pair_idx[1,]], g[pair_idx[2,]]])
    })
    
    k <- which.max(score)
    return(list(stations = st[cmb[, k]], score = score[k], D = D))
  }
  
  # --- greedy heuristic -----------------------------------------------------
  # pick the first two
  if (length(fixed_idx) == 0) {
    i1 <- which.max(rowSums(D))
    i2 <- which.max(D[i1, ])
  } else if (length(fixed_idx) == 1) {
    i1 <- fixed_idx
    i2 <- which.max(D[i1, ])
  } else if (length(fixed_idx) >= 2) {
    i1 <- fixed_idx[1]
    i2 <- fixed_idx[2]
  }
  
  chosen <- c(i1, i2)
  
  # pick third
  min_dist_to_chosen <- pmin(D[, i1], D[, i2])
  min_dist_to_chosen[chosen] <- -Inf
  if (length(fixed_idx) >= 3) {
    i3 <- fixed_idx[3]
  } else {
    i3 <- which.max(min_dist_to_chosen)
  }
  chosen <- c(chosen, i3)
  
  # pick fourth: maximize min distance to already chosen
  tmp <- apply(D[, chosen, drop = FALSE], 1, min)
  tmp[chosen] <- -Inf
  i4 <- which.max(tmp)
  chosen <- c(chosen, i4)
  
  # score = sum of all 6 distances
  idx <- combn(chosen, 2)
  score <- sum(D[idx[1,], idx[2,]])
  
  list(stations = st[chosen], score = score, D = D)
}

# ------------------------------
#### Data Exploration ####
# ------------------------------
Q <- station_quantile_matrix(df, var = "wind_speed_of_gust", by = "location")

# Plot 4 most heterogeneous stations
sel <- top4_most_heterogeneous(Q, exact = TRUE)
sel$stations  # <- the 4 most heterogeneous stations
s4 <- forecasts_wind_gust
s4 <- df[df$location %in% sel$stations, , drop = FALSE]
ggplot(s4, aes(x = wind_speed_of_gust, fill = location)) +
  geom_density(alpha = 0.35, adjust = 1) +
  labs(title = "Windgust Distribution of Selected Heterogeneous Stations",
       x = "Wind speed of gust", y = "Density") +
  theme_minimal()


# Plot 4 most heterogeneous stations including the one next to Karlsruhe and Sylt (northern most station)
sel_ka <- top4_most_heterogeneous(Q, exact = TRUE, fixed_stations=c(10731, 10020))
sel_ka$stations  # <- the 4 most heterogeneous stations (Sylt 10020, Brocken 10453, Rheinstetten 10731, Garmisch-Partenkirchen 10963)

station_names <- c(
  "10020" = "Sylt",
  "10453" = "Brocken",
  "10731" = "Rheinstetten",
  "10963" = "Garmisch-Partenkirchen"
)

s4_ka <- df %>%
  filter(location %in% sel_ka$stations) %>%
  mutate(station_name = station_names[as.character(location)])

ggplot(s4_ka, aes(x = wind_speed_of_gust, fill = station_name)) +
  geom_density(alpha = 0.35, adjust = 4) +
  labs(
    title = "Wind Gust Distribution of Selected Stations",
    x = "Wind speed of gust", y = "Density",
    fill = "Station"
  ) +
  theme_minimal(base_family = "Franklin") +
  theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12)
  )

# 10020 Sylt (Germany’s northernmost island, directly facing the North Sea):
# The area is extremely exposed to Atlantic westerlies, with almost no topographic shielding.
# Wind speeds are consistently high, and gusts pick up quickly even without storms.
# In the density plot:
#  The pink distribution starts low (few calm periods).
#  It has a broad main peak around moderate gust speeds (6–10 m/s).
#  The long right tail shows that very strong gusts (20–30+ m/s) are common.

# 10453 Brocken (mountain summit station 3with no surface roughness):
# Highest peak in northern Germany and one of the windiest places in the entire country.
# Completely exposed at summit level.
# Known for severe storms, hurricane-force winds, and frequent low-altitude jet interactions.
# In the density plot:
#   The yellow distribution has:
#   very wide spread (gust speeds vary enormously),
#   the heaviest tail — Brocken frequently produces the strongest gusts of all four,
#   a broad plateau rather than a sharp peak → gusts are variable and turbulent.

# 10731 Karlsruhe-Rheinstetten (Upper Rhine Valley / lowland, inland)
# Located in the Upper Rhine Valley in southwestern Germany.
# Sheltered by the Black Forest to the east and Vosges Mountains to the west.
# Known for low to moderate winds; strong gusts occur mostly during deep cyclones.
# In the density plot:
#   The blue density is centered lower (~7–10 m/s).
#   It has a much shorter tail → extreme gusts >20 m/s are rare.
#   The distribution is narrower → less variation than coastal or mountain regions.

# 10963 – Garmisch-Partenkirchen (sheltered Alpine valley station)
#   Surrounding mountains block much of the synoptic wind.
#   Gusts tend to be terrain-driven (katabatic winds, valley channeling) and not extreme.
#   In the density plot:
#   The purple curve is left-shifted → mostly weaker gusts (3–9 m/s).
#   It has very little right tail → strong storms rarely reach the valley floor.
#   The peak is sharp → wind conditions are stable and predictable.

# Plot 4 most heterogeneous stations emphasizing the importance of the tails of the distribution
# baseline weight 1, heavier for probs >= 0.9
probs <- seq(0.01, 0.99, by = 0.01)
weights <- ifelse(probs >= 0.9, 10, 1)

sel_tail <- top4_most_heterogeneous(
  Q,
  exact = TRUE,
#  fixed_stations=c(10731),  # include Rheinstetten
  weights = weights
)

sel_tail$stations

s4_tail <- df[df$location %in% sel_tail$stations, , drop = FALSE]
ggplot(s4_tail, aes(x = wind_speed_of_gust, fill = location)) +
  geom_density(alpha = 0.35, adjust = 1) +
  labs(title = "Windgust Distribution of Selected Heterogeneous Stations including Rheinstetten",
       x = "Wind speed of gust", y = "Density") +
  theme_minimal()


# Plot 4 random stations
set.seed(123)  # for reproducibility
random_stations <- sample(setdiff(unique(df$location), sel$stations), 4)
# subset the data
s4_random <- df[df$location %in% random_stations, , drop = FALSE]
ggplot(s4_random, aes(x = wind_speed_of_gust, fill = location)) +
  geom_density(alpha = 0.35, adjust = 1) +
  labs(title = "Windgust Distribution of Randomly Selected Stations",
       x = "Wind speed of gust", y = "Density") +
  theme_minimal()


#### Plot wind time series of 4 selected stations ####

ggplot(s4_ka, aes(x = obs_tm, y = wind_speed_of_gust)) +
  geom_line(aes(color = station_name), alpha = 0.6, linewidth = 0.4)+
  facet_wrap(~ station_name, ncol = 1, scales = "fixed") + 
  labs(
    title = "Wind Gusts in m/s Over Time for Selected Stations",
    x = "Year",
    y = "Wind speed of gust"
  ) +
  theme_minimal(base_family = "Franklin") +  
  theme(
    strip.text = element_text(size = 25),
    plot.title = element_text(size = 25),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 20),
    legend.position = "none"   # None needed—station names are in facet labels
  )


# ------------------------------
#### Get summary statistics ####
# ------------------------------
stations_to_check <- c("10020", "10453", "10731", "10963")

summary_obs_long <- df %>%
  filter(location %in% stations_to_check) %>%
  group_by(location) %>%
  summarise(across(
    all_of(obs_vars),
    list(
      missing = ~ sum(is.na(.))/length(.)*100,
      min     = ~ min(., na.rm = TRUE),
      max     = ~ max(., na.rm = TRUE)
    )
  )) %>%
  pivot_longer(
    cols = -location,
    names_to      = c("variable", "stat"),
    names_pattern = "(.+)_([^_]+)$",  # everything_before_last_ , last_piece
    values_to     = "value"
  )

print(n = 100, summary_obs_long)

summary_obs_wide <- summary_obs_long %>%
  mutate(
    stat = factor(stat, levels = c("min", "max", "missing"))
  ) %>%
  arrange(variable, stat) %>%           # group mins / maxes / missing together
  pivot_wider(
    id_cols  = c(variable, stat),       # each row = variable × stat
    names_from  = location,             # columns = stations
    values_from = value
  )

summary_obs_wide

