library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

source("~/scripts/QRN_EQRN.R")
source("~/scripts/iid_RF_EQRN.R")

cov_dir <- "~/coverage"

station_names <- c(
  "10020" = "Sylt",
  "10453" = "Brocken",
  "10731" = "Rheinstetten",
  "10963" = "Garmisch-Partenkirchen"
)

# ---------- helper to load a single .RData file and extract 1 numeric object ----------
read_one_value <- function(path, obj_name) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (!exists(obj_name, envir = e)) stop("No object '", obj_name, "' in ", path)
  as.numeric(get(obj_name, envir = e))
}

# ---------- generic reader: parse setup + station from filename ----------
parse_setup_station <- function(fname) {
  tibble(
    setup   = str_match(fname, "^(iid|rnn)_")[,2],
    station = str_match(fname, "_([0-9]+)$")[,2]
  )
}

# ---------- 1) EQRN coverage ----------
files_cov_eqrn <- list.files(cov_dir, pattern = "^(iid|rnn)_coverage_[0-9]+$", full.names = TRUE)

cov_eqrn_long <- map_dfr(files_cov_eqrn, function(path) {
  fname <- basename(path)
  meta <- parse_setup_station(fname)
  tibble(
    station  = meta$station,
    setup    = meta$setup,
    tau      = 0.99,
    coverage = read_one_value(path, "coverage")
  )
})

# ---------- 2) EQRN quantile score ----------
files_qs_eqrn <- list.files(cov_dir, pattern = "^(iid|rnn)_qs_eqrn_[0-9]+$", full.names = TRUE)

qs_eqrn_long <- map_dfr(files_qs_eqrn, function(path) {
  fname <- basename(path)
  meta <- parse_setup_station(fname)
  tibble(
    station = meta$station,
    setup   = meta$setup,
    tau     = 0.99,
    qs      = read_one_value(path, "qs_eqrn")  # <- change if you saved under a different object name
  )
})

# ---------- 3) QRN coverage ----------
files_cov_qrn <- list.files(cov_dir, pattern = "^(iid|rnn)_coverage_qrn_[0-9]+$", full.names = TRUE)

cov_qrn_long <- map_dfr(files_cov_qrn, function(path) {
  fname <- basename(path)
  meta <- parse_setup_station(fname)
  tibble(
    station  = meta$station,
    setup    = meta$setup,
    tau      = 0.99,
    coverage = read_one_value(path, "coverage_qrn")  # <- change if object name differs
  )
})

# ---------- 4) QRN quantile score ----------
files_qs_qrn <- list.files(cov_dir, pattern = "^(iid|rnn)_qs_qrn_[0-9]+$", full.names = TRUE)

qs_qrn_long <- map_dfr(files_qs_qrn, function(path) {
  fname <- basename(path)
  meta <- parse_setup_station(fname)
  tibble(
    station = meta$station,
    setup   = meta$setup,
    tau     = 0.99,
    qs      = read_one_value(path, "qs_qrn")   # <- change if you saved under a different object name
  )
})

# ---------- Build the two final tables ----------
make_wide_table <- function(cov_long, qs_long) {
  full_join(cov_long, qs_long, by = c("station", "setup", "tau")) %>%
    mutate(station_name = station_names[station]) %>%
    select(station, station_name, setup, tau, coverage, qs) %>%
    pivot_wider(
      names_from = setup,
      values_from = c(coverage, qs),
      names_sep = "_"
    ) %>%
    arrange(station)
}

table_eqrn <- make_wide_table(cov_eqrn_long, qs_eqrn_long)
table_qrn  <- make_wide_table(cov_qrn_long,  qs_qrn_long)

table_eqrn
table_qrn
