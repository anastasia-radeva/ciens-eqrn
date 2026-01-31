# CIENS DATA PROCESSING

## IMPORTANT NOTE: You require at least the folder "run0" and "run12" and the folder "observations" for this code to work. 

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(purrr)
library(forcats)
library(matrixStats)


# ------------------------------------------------------------------------------
# Initialize
# ------------------------------------------------------------------------------

# Path to Github repository functions
git_path <- "~/"
# Set working directory
setwd(git_path)
# Source R file for initializing
source(file = paste0(getwd(), "/init_file.R"))  # Specify data path here

# ------------------------------------------------------------------------------
# Construct data.frame from netcdf files
# ------------------------------------------------------------------------------

# Isolate selected stations
loc_vec <- loc_data$station_id[is.element(loc_data$station_id, c("10731", "10020", "10963", "10453"))]
# Restrict to initialization times at 12  and 00 UTC for forecasts
tm_vec <- init_vec[hour(init_vec) == 12 | hour(init_vec) == 00]
tm_vec <- tm_vec[year(tm_vec) >= 2010 & year(tm_vec) <= 2023]
# Wind gust forecasts variable
met_vars <- c("VMAX_10M")

# Generate data.frame df for all initialization times and all forecast horizons

file_name <- file.path("~/forecasts/forecasts_wind_gust")

if (file.exists(file_name)) {
  load(file = "~/forecasts/forecasts_wind_gust")

} else {
  forecasts_wind_gust <- bind_rows(lapply(tm_vec, function(x) get_init(tm = x,
                                                       dir_path = data_path,
                                                       met_vars = met_vars,
                                                       location_vec = loc_vec,
                                                       step_vec = c(0:21),
                                                       ens_vec = ens_vec)))
  
  save(forecasts_wind_gust, file = "~/forecasts/forecasts_wind_gust")
}

# ------------------------------------------------------------------------------
# Construct separate data frames for each station, 
# compute ensemble summary statistics
# ------------------------------------------------------------------------------

# Helper function
add_ensemble_summary <- function(data,
                                 station_id,
                                 ens_prefix = "^VMAX_10M_",
                                 obs_var = "wind_speed_of_gust",
                                 keep_ens_cols = TRUE) {
  
  # Filter station and remove missing observations
  station_df <- data %>%
    filter(location == station_id, !is.na(.data[[obs_var]]))
  
  # Identify ensemble member columns
  ens_cols <- grep(ens_prefix, names(station_df), value = TRUE)
  if (length(ens_cols) == 0) {
    stop("No ensemble columns found for prefix: ", ens_prefix)
  }
  
  # Convert selected ensemble columns to matrix (needed for matrixStats speed)
  ens_mat <- as.matrix(station_df[, ens_cols, drop = FALSE])
  
  # Add row-wise mean + sd
  station_df <- station_df %>%
    mutate(
      ens_mean = rowMeans(ens_mat, na.rm = TRUE),
      ens_sd   = matrixStats::rowSds(ens_mat, na.rm = TRUE)
    )
  
  # Optional: drop raw ensemble columns to save memory
  if (!keep_ens_cols) {
    station_df <- station_df %>% select(-all_of(ens_cols))
  }
  
  return(station_df)
}


# Apply to all four stations
stations <- c("10020", "10453", "10731", "10963")

station_dfs <- list()

for (st in stations) {
  station_dfs[[st]] <- add_ensemble_summary(
    data = forecasts_wind_gust,
    station_id = st,
    ens_prefix = "^VMAX_10M_",
    obs_var = "wind_speed_of_gust",
    keep_ens_cols = TRUE
  )
}

forecasts_10020 <- station_dfs[["10020"]]
forecasts_10453 <- station_dfs[["10453"]]
forecasts_10731 <- station_dfs[["10731"]]
forecasts_10963 <- station_dfs[["10963"]]
