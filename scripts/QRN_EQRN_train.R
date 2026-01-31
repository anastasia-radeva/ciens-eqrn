# QRN-EQRN Model Training
# Example script for the Brocken station.
# If you want to train for another station, 
# replace "10453" with the identifier of your desired station.

library(dplyr)
library(EQRN)
library(grf)
library(torch)
source("~/scripts/patches/GPD_excess_probability_patch.R")
source("~/scripts/patches/QRN_patch.R")
source("~/scripts/patches/EQRN_seq_patch.R")
source("~/scripts/data_preprocessing.R")

source(file = paste0(getwd(), "/init_file.R"))

# ------------------------------------------------------------
# Build 12-step forecast sequence mapped to obs at 00/12 UTC
# ------------------------------------------------------------
forecasts_10453 <- forecasts_10453 %>% filter(step %in% 1:12)

# ------------------------------------------------------------
# Train/test split
# ------------------------------------------------------------

cutoff_date <- as.POSIXct("2021-06-01 01:00:00", tz = "UTC")
train_df <- forecasts_10453 %>% filter(obs_tm < cutoff_date)
test_df  <- forecasts_10453 %>% filter(obs_tm >= cutoff_date)

rownames(train_df) <- train_df$obs_tm
rownames(test_df) <- test_df$obs_tm

# Duplicate target values across all steps of series
k <- 12
# train_df <- train_df %>%
#   arrange(location, init_tm, step) %>%
#   group_by(location, init_tm) %>%
#   mutate(wind_speed_of_gust = wind_speed_of_gust[step == k][1]) %>%
#   ungroup()
# 
# test_df <- test_df %>%
#   arrange(location, init_tm, step) %>%
#   group_by(location, init_tm) %>%
#   mutate(wind_speed_of_gust = wind_speed_of_gust[step == k][1]) %>%
#   ungroup()

x_vars <- c("obs_tm", "step", "ens_mean", "ens_sd")
y_vars <- c("obs_tm", "step", "wind_speed_of_gust")

# Remove missing
train_df <- train_df %>% filter(if_all(all_of(c(x_vars, y_vars)), ~ is.finite(.)))
test_df  <- test_df  %>% filter(if_all(all_of(c(x_vars, y_vars)), ~ is.finite(.)))

train_df <- train_df %>%
  group_by(location, init_tm) %>%
  filter(setequal(unique(step), 1:k)) %>%
  ungroup()

test_df <- test_df %>%
  group_by(location, init_tm) %>%
  filter(setequal(unique(step), 1:k)) %>%
  ungroup()

X_train <- as.matrix(train_df[, x_vars, drop = FALSE])
y_train <- as.matrix((train_df)[, y_vars, drop = FALSE])

X_test  <- as.matrix(test_df[, c("ens_mean","ens_sd"), drop = FALSE])
y_test  <- as.matrix((test_df)[, "wind_speed_of_gust", drop = FALSE])


# Validation split
cutoff_date_val <- as.POSIXct("2019-02-22 00:00:00", tz = "UTC")  # around 80% of training data

idx_tr  <- X_train[, "obs_tm"] < cutoff_date_val

X_tr  <- X_train[idx_tr, c("ens_mean", "ens_sd"), drop = FALSE]
X_val <- X_train[!idx_tr, c("ens_mean", "ens_sd"), drop = FALSE]
X_tr  <- apply(X_tr,  2, as.numeric)
X_val <- apply(X_val, 2, as.numeric)

y_tr   <- as.numeric(y_train[idx_tr, "wind_speed_of_gust", drop = FALSE])
y_val  <- as.numeric(y_train[!idx_tr, "wind_speed_of_gust", drop = FALSE])


# ------------------------------------------------------------
# Intermediate quantile model (recurrent QRN)
# ------------------------------------------------------------
interm_lvl <- 0.8
levels_predict <- 0.99

model_dir  <- "~/fits"
qrn_model_name <- "fit_QRN_12stride_10453"

qrn_model_file <- file.path(model_dir, qrn_model_name)

if (file.exists(qrn_model_file)) {
  fit_QRN_ts <- EQRN_load(path = model_dir, name = qrn_model_name, device = default_device())
  fit_QRN_ts$fit_nn <- fit_QRN_ts$fit_nn$to(device = default_device())
  fit_QRN_ts$fit_nn$eval()
} else {
  fit_QRN_ts <- QRN_seq_fit(
    X=X_tr, Y=y_tr, q_level=interm_lvl,
    hidden_size=32, num_layers=2, rnn_type="lstm",
    p_drop=0.1,
    learning_rate=3e-4, L2_pen=1e-6,
    seq_len=12, stride = 12, scale_features=TRUE,
    n_epochs=200, batch_size=512,
    X_valid=X_val, Y_valid=y_val,
    lr_decay=0.5, patience_decay=5, min_lr=1e-5, patience_stop=20,
    tol=1e-4, fold_separation=NULL, warm_start_path=NULL, patience_lag=5, optim_met="adam",
    seed=2, verbose=2, device=default_device()
  )
  
  EQRN_save(fit_QRN_ts, path = model_dir, name = qrn_model_name, no_warning = TRUE)
}


# ------------------------------------------------------------
# Tail model (recurrent EQRN)
# ------------------------------------------------------------

# Construct out-of-bag intermediate quantiles on the training set.
intermediateq_train <- QRN_seq_predict(fit_QRN_ts, X_tr, y_tr, q_level=interm_lvl, crop_predictions=FALSE, device=default_device())

# Fit the tail model -------------------------------------------------------
quant_val <- QRN_seq_predict(
  fit_QRN_ts, X_val, y_val,
  q_level = interm_lvl,
  crop_predictions = FALSE
)


eqrn_model_name <- "fit_EQRN_12stride_10453"

eqrn_model_file <- file.path(model_dir, eqrn_model_name)

if (file.exists(eqrn_model_file)) {
  fit_EQRN_ts <- EQRN_load(path = model_dir, name = eqrn_model_name, device = default_device())
  fit_EQRN_ts$fit_nn <- fit_EQRN_ts$fit_nn$to(device = default_device())
  fit_EQRN_ts$fit_nn$eval()
} else {
  fit_EQRN_ts <- EQRN_fit_seq(
    X_tr, y_tr,
    intermediate_quantiles=intermediateq_train,
    interm_lvl=interm_lvl,
    shape_fixed=FALSE,
    hidden_size=64, num_layers=1, rnn_type="lstm",
    p_drop=0.1,
    intermediate_q_feature=FALSE,
    learning_rate=1e-4, L2_pen=1e-6,
    seq_len=12, stride = 12, 
    shape_penalty=0,
    scale_features=TRUE,
    n_epochs=300, batch_size=512,
    X_valid=X_val, y_valid=y_val, quant_valid=quant_val,
    lr_decay=0.5, patience_decay=5, min_lr=1e-6, patience_stop=30,
    tol=1e-5,
    orthogonal_gpd=TRUE,
    patience_lag=1, fold_separation=NULL, optim_met="adam",
    seed=2, verbose=2, device=default_device()
  )
  
  EQRN_save(fit_EQRN_ts, path = model_dir, name = eqrn_model_name, no_warning = TRUE)
}


# ------------------------------------------------------------
# Predict extreme quantile
# ------------------------------------------------------------
# Predict intermediate test quantiles using the intermediate model.
intermediateq_test_ts <- QRN_seq_predict(fit_QRN_ts, X_test, y_test, q_level=interm_lvl, crop_predictions=FALSE, device=default_device())

qpred_eqrn_ts <- EQRN_predict_seq(
  fit_EQRN_ts,
  X = X_test,
  Y = y_test,
  prob_lvls_predict = levels_predict,
  intermediate_quantiles = intermediateq_test_ts,
  interm_lvl = interm_lvl
)

# ------------------------------------------------------------
# Exceedance probability at empirical 0.99 test quantile
# ------------------------------------------------------------
emp_quantile <- as.numeric(quantile(y_test, 0.99, na.rm = TRUE))

ppred_eqrn_ts <- EQRN_excess_probability_seq(
  val = emp_quantile,
  fit_eqrn = fit_EQRN_ts,
  X = X_test,
  Y = y_test,
  intermediate_quantiles = intermediateq_test_ts,
  interm_lvl = interm_lvl,
  proba_type = "excess"
)

# Check
i <- which.max(y_test)
i
ppred_eqrn_ts[i]

# ------------------------------------------------------------
# Compute baseline using quantile regression network
# ------------------------------------------------------------
extreme_qrn_model_name <- "fit_QRN_extreme_12stride_10453"

extreme_qrn_model_file <- file.path(model_dir, extreme_qrn_model_name)

if (file.exists(extreme_qrn_model_file)) {
  fit_QRN_ts_extreme <- EQRN_load(path = model_dir, name = extreme_qrn_model_name, device = default_device())
  fit_QRN_ts_extreme$fit_nn <- fit_QRN_ts_extreme$fit_nn$to(device = default_device())
  fit_QRN_ts_extreme$fit_nn$eval()
} else {
  fit_QRN_ts_extreme <- QRN_seq_fit(
    X=X_tr, Y=y_tr, q_level=levels_predict,
    hidden_size=32, num_layers=2, rnn_type="lstm",
    p_drop=0.1,
    learning_rate=3e-4, L2_pen=1e-6,
    seq_len=12, stride = 12, scale_features=TRUE,
    n_epochs=200, batch_size=512,
    X_valid=X_val, Y_valid=y_val,
    lr_decay=0.5, patience_decay=5, min_lr=1e-5, patience_stop=20,
    tol=1e-4, fold_separation=NULL, warm_start_path=NULL, patience_lag=5, optim_met="adam",
    seed=2, verbose=2, device=default_device()
  )
  
  EQRN_save(fit_QRN_ts_extreme, path = model_dir, name = extreme_qrn_model_name, no_warning = TRUE)
}

# Compute quantiles
qpred_qrn_ts <- QRN_seq_predict(fit_QRN_ts_extreme, X_test, y_test, q_level=levels_predict, crop_predictions=FALSE, device=default_device())
