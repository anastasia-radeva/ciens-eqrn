# EVALUATION PLOTS

# Example script for the Brocken station.
# If you want to evaluate for another station, 
# replace "10453" with the identifier of your desired station.

# -----------------------------
# Dependencies
# -----------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(patchwork)

source("~/scripts/QRN_EQRN_train.R")

Sys.setlocale("LC_TIME", "en_US.UTF-8")  # set plot language to English

# Helper filter function
is_00_12 <- function(x) {
  hr <- lubridate::hour(x)
  hr %in% c(0, 12)
}

# -----------------------------
# Ensemble member columns
# -----------------------------
ens_cols <- grep("^VMAX_10M_", names(forecasts_10453), value = TRUE)
if (length(ens_cols) == 0) stop("No ensemble columns found matching '^VMAX_10M_'")

# -----------------------------
# Prepare obs + predicted quantile vectors
# -----------------------------
y_vec <- as.numeric(y_test)

if (is.matrix(qpred_eqrn_ts)) {
  q_vec <- as.numeric(qpred_eqrn_ts[, 1])
} else {
  q_vec <- as.numeric(qpred_eqrn_ts)
}

if (is.matrix(qpred_qrn_ts)) {
  q_qrn_vec <- as.numeric(qpred_qrn_ts[, 1])
} else {
  q_qrn_vec <- as.numeric(qpred_qrn_ts)
}

stopifnot(nrow(test_df) == length(y_vec))
stopifnot(nrow(test_df) == length(q_vec))


# -----------------------------
# Build table for plotting
# -----------------------------
ts_df <- test_df %>%
  transmute(
    obs_tm = obs_tm,
    y_obs  = y_vec,
    q_ext  = q_vec,
    q_qrn  = q_qrn_vec
  ) %>%
  filter(is_00_12(obs_tm))

# -----------------------------
# Global coverage on test set
# -----------------------------
coverage_qrn <- mean(ts_df$y_obs <= ts_df$q_qrn, na.rm = TRUE)
save(coverage_qrn, file = "~/coverage/rnn_coverage_qrn_10453")

coverage <- mean(ts_df$y_obs <= ts_df$q_ext, na.rm = TRUE)
save(coverage, file = "~/coverage/rnn_coverage_10453")
n_exceed <- sum(ts_df$y_obs > ts_df$q_ext, na.rm = TRUE)
n_total <- sum(is.finite(ts_df$y_obs) & is.finite(ts_df$q_ext))

coverage_txt <- paste0(
  "Global Coverage = ", sprintf("%.4f", coverage), "\n",
  "Exceedances = ", n_exceed, " / ", n_total
)

# -----------------------------
# Quantile score
# -----------------------------
quantile_score <- function(y, q, tau) {
  y <- as.numeric(y)
  q <- as.numeric(q)
  ok <- is.finite(y) & is.finite(q)
  y <- y[ok]; q <- q[ok]
  mean((tau - (y <= q)) * (y - q))
}

qs_qrn  <- quantile_score(ts_df$y_obs, ts_df$q_qrn, levels_predict)
qs_eqrn <- quantile_score(ts_df$y_obs, ts_df$q_ext, levels_predict)

save(qs_qrn,  file = "~/coverage/rnn_qs_qrn_10453")
save(qs_eqrn, file = "~/coverage/rnn_qs_eqrn_10453")

# -----------------------------
# Compute peak
# -----------------------------
# Identify peak + window
i_peak <- which.max(ts_df$y_obs)
t0 <- ts_df$obs_tm[i_peak]

window_hours <- 240
t_min <- t0 - hours(window_hours)
t_max <- t0 + hours(window_hours)

ts_win <- ts_df %>%
  filter(obs_tm >= t_min, obs_tm <= t_max)

# -----------------------------
# Compute ensemble MAX per obs_tm within the same window
# -----------------------------
ens_max_win <- forecasts_10453 %>%
  filter(obs_tm >= t_min, obs_tm <= t_max) %>%
  filter(is_00_12(obs_tm)) %>%
  select(obs_tm, all_of(ens_cols)) %>%
  mutate(ens_max = apply(across(all_of(ens_cols)), 1, max, na.rm = TRUE)) %>%
  select(obs_tm, ens_max)

# -----------------------------
# Long format for observation + EQRN quantile
# -----------------------------
lines_long <- ts_win %>%
  select(obs_tm, y_obs, q_ext) %>%
  pivot_longer(cols = c(y_obs, q_ext),
               names_to = "series",
               values_to = "value") %>%
  mutate(series = recode(series,
                         y_obs = "Observation",
                         q_ext = paste0("EQRN ", levels_predict, "-Quantile")))

# -----------------------------
# Plot
# -----------------------------

# Calibration Plot -----------------------------
eqrn_line <- lines_long %>% filter(series != "Observation")
obs_line  <- ts_win %>% select(obs_tm, y_obs)

legend_levels <- c(
  "Observation",
  "Empirical 0.99-Quantile",
  "Forecast Ensemble Max",
  "EQRN 0.99-Quantile"
)

legend_colors <- c(
  "Observation" = "#4E79A7",
  "Empirical 0.99-Quantile" = "#D55E00",
  "Forecast Ensemble Max" = "#59A14F",
  "EQRN 0.99-Quantile" = "#7A0019"
)

color_scale <- scale_color_manual(
  values = legend_colors,
  breaks = legend_levels,
  drop = FALSE
)

legend_guide <- guides(color = guide_legend(title = NULL))

base_theme <- theme_minimal(base_size = 13, base_family = "Franklin Gothic") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

y_limits <- coord_cartesian(ylim = c(0, 70))

# Just observations
p1 <- ggplot() +
  geom_line(
    data = obs_line,
    aes(x = obs_tm, y = y_obs, color = "Observation"),
    linewidth = 1.1
  ) +
  labs(
    title = "Wind Gusts in Brocken",
    subtitle = paste("Window centered on peak gust at", format(t0, "%Y-%m-%d %H:%M:%S")),
    x = "Observation time",
    y = "Wind gust"
  ) +
  color_scale +
  legend_guide +
  base_theme +
  y_limits

# Observations and empirical quantile
p2 <- p1 +
  geom_hline(
    aes(yintercept = emp_quantile, color = "Empirical 0.99-Quantile"),
    linetype = "dotted",
    linewidth = 1.1
  )

# Observations, empirical quantile and ensemble max
p3 <- p2 +
  geom_line(
    data = ens_max_win,
    aes(x = obs_tm, y = ens_max, color = "Forecast Ensemble Max"),
    linewidth = 1.1,
    alpha = 0.8
  )


# Observations, empirical quantile, ensemble max and EQRN Quantile
p4 <- p3 +
  geom_line(
    data = eqrn_line,
    aes(x = obs_tm, y = value, color = "EQRN 0.99-Quantile"),
    linewidth = 1.1
  )

p5 <- p1 +
  geom_line(
    data = eqrn_line,
    aes(x = obs_tm, y = value, color = "EQRN 0.99-Quantile"),
    linewidth = 1.1
  ) +
  # annotate coverage 
  annotate( "label", x = min(ts_win$obs_tm), 
            y = max(c(ts_win$y_obs, 
                      ens_max_win$ens_max), 
                    na.rm = TRUE), 
            label = coverage_txt, 
            hjust = 0, vjust = 1, 
            size = 3.8, label.size = 0.25, 
            fill = "white", alpha = 0.9 )

# Display plots
p1
p2
p3
p4
p5

# Exceedance Probability Plot -----------------------------
# --- build full df first (aligned vectors)
uncond_p_exc <- mean(y_test > emp_quantile, na.rm = TRUE) 

plot_df <- test_df %>%
  transmute(
    obs_tm = obs_tm,
    y_obs  = y_vec,
    p_exc  = ppred_eqrn_ts,
    exc_ratio = p_exc / uncond_p_exc
  ) %>%
  filter(is_00_12(obs_tm))

plot_win <- plot_df %>%
  filter(obs_tm >= t_min, obs_tm <= t_max)

# (optional) horizontal threshold line if you want it
# emp_quantile <- quantile(plot_df$y_obs, 0.99)


# --- Top panel: observations
p_top <- ggplot(plot_win, aes(x = obs_tm)) +
  geom_line(aes(y = y_obs, color = "Observation"), linewidth = 1) +
  geom_hline(yintercept = emp_quantile, linetype = "dashed", color = "grey30") +
  geom_vline(xintercept = t0, linetype = "dashed", color = "black") +
  labs(
    y = "Wind gust",
    x = "",
    title = "Predicted Exceedance Probability Ratio in Brocken",
    subtitle = paste("Window centered on peak gust at", format(t0, "%Y-%m-%d %H:%M:%S"))
  ) +
  scale_color_manual(values = c(
    "Observation" = "#4E79A7"
  )) +
  guides(color = guide_legend(title = NULL)) +
  base_theme +
  coord_cartesian(xlim = c(t_min, t_max)) +
  y_limits

# --- Bottom panel: exceedance probability
p_bottom <- ggplot(plot_win, aes(x = obs_tm)) +
  geom_line(aes(y = exc_ratio, color = "P(Y > Empirical 0.99-Quantile | X) / P(Y > Empirical 0.99-Quantile)"), linewidth = 1) +
  geom_vline(xintercept = t0, linetype = "dashed", color = "black") +
  labs(
    x = "Observation time",
    y = "Probability ratio"
  ) +
  scale_color_manual(values = c(
    "P(Y > Empirical 0.99-Quantile | X) / P(Y > Empirical 0.99-Quantile)" = "#8B0A50"
  )) +
  guides(color = guide_legend(title = NULL)) +
  base_theme +
  coord_cartesian(xlim = c(t_min, t_max)) +
  ylim(0, 25)


# --- Combine
p_top / p_bottom + plot_layout(heights = c(2, 1))

