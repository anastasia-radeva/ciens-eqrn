# GRF-EQRN Model Training
# Example script for the Brocken station.
# If you want to train for another station, 
# replace "10453" with the identifier of your desired station.

library(dplyr)
library(lubridate)
library(EQRN)
library(grf)
library(torch)
source("~/scripts/patches/GPD_excess_probability_patch.R")
source("~/scripts/data_preprocessing.R")
source(file = paste0(getwd(), "/init_file.R"))

# -----------------------------
# PARAMETERS
# -----------------------------
step_target  <- 12
cutoff_date  <- as.POSIXct("2021-06-01", tz = "UTC")

x_vars <- c("obs_tm","step","ens_mean", "ens_sd")
y_var  <- "wind_speed_of_gust"

interm_lvl      <- 0.8
levels_predict  <- 0.99

# -----------------------------
# FILTER: init in {00,12}, step=12, obs in {00,12}
# -----------------------------
df_step12 <- forecasts_10453 %>%
  filter(
    step == step_target,
    hour(obs_tm)  %in% c(0, 12)
  )

# -----------------------------
# TRAIN / TEST split by obs time
# -----------------------------
train_df <- df_step12 %>% filter(obs_tm < cutoff_date)
test_df  <- df_step12 %>% filter(obs_tm >= cutoff_date)

# remove missing in X or Y
train_df <- train_df %>%
  filter(if_all(all_of(c(x_vars, y_var)), ~ is.finite(.)))

test_df <- test_df %>%
  filter(if_all(all_of(c(x_vars, y_var)), ~ is.finite(.)))

# -----------------------------
# MATRICES
# -----------------------------
X_train <- as.matrix(train_df[, x_vars, drop = FALSE])
y_train <- as.matrix(train_df[, y_var,  drop = FALSE])

X_test  <- as.matrix(test_df[,  x_vars, drop = FALSE])
y_test  <- as.matrix(test_df[,  y_var,  drop = FALSE])

cat("Train rows:", nrow(X_train), "\n")
cat("Test rows :", nrow(X_test),  "\n")

# -----------------------------
# VALIDATION SPLIT (IID)
# -----------------------------
set.seed(2)

n   <- nrow(X_train)
idx <- sample.int(n)

cut <- floor(0.8 * n)
train_idx <- idx[1:cut]
val_idx   <- idx[(cut+1):n]

X_tr  <- X_train[train_idx, , drop = FALSE]
y_tr  <- y_train[train_idx, , drop = FALSE]

X_val <- X_train[val_idx, , drop = FALSE]
y_val <- y_train[val_idx, , drop = FALSE]

# -----------------------------
# STEP 1: INTERMEDIATE MODEL (GRF)
# -----------------------------
fit_grf <- grf::quantile_forest(X_tr, y_tr, num.trees = 1000, seed = 2)

intermediateq_train <- predict(fit_grf, newdata = NULL, quantiles = interm_lvl)$predictions
quant_val           <- predict(fit_grf, newdata = X_val, quantiles = interm_lvl)$predictions
intermediateq_test  <- predict(fit_grf, newdata = X_test, quantiles = interm_lvl)$predictions

# -----------------------------
# STEP 2: FIT IID EQRN
# -----------------------------
model_dir  <- "~/fits"
eqrn_model_name <- "fit_EQRN_iid_10453"

eqrn_model_file <- file.path(model_dir, eqrn_model_name)

if (file.exists(eqrn_model_file)) {
  fit_EQRN <- EQRN_load(path = model_dir, name = eqrn_model_name, device = default_device())
  fit_EQRN$fit_nn <- fit_EQRN$fit_nn$to(device = default_device())
  fit_EQRN$fit_nn$eval()
} else {
  fit_EQRN <- EQRN_fit(
    X = X_tr,
    y = y_tr,
    intermediate_quantiles = intermediateq_train,
    interm_lvl = interm_lvl,
    intermediate_q_feature = FALSE,
    shape_fixed = FALSE,
    
    # Good efficiency/performance tradeoff for only 2 features:
    net_structure = c(32, 16),
    hidden_fct = torch::nnf_relu,      # faster & typically better than sigmoid
    p_drop = 0.1,
    
    learning_rate = 1e-3,
    L2_pen = 1e-6,
    n_epochs = 250,
    batch_size = 512,
    
    X_valid = X_val,
    y_valid = y_val,
    quant_valid = quant_val,
    
    lr_decay = 0.5,
    patience_decay = 10,
    min_lr = 1e-6,
    patience_stop = 30,
    tol = 1e-5,
    
    orthogonal_gpd = TRUE,
    seed = 2,
    verbose = 2,
    device = default_device()
  )
  
  EQRN_save(fit_EQRN, path = model_dir, name = eqrn_model_name, no_warning = TRUE)
}

# -----------------------------
# STEP 3: EXTREME QUANTILES
# -----------------------------
qpred_eqrn <- EQRN_predict(
  fit_EQRN,
  X = X_test,
  prob_lvls_predict = levels_predict,
  intermediate_quantiles = intermediateq_test,
  interm_lvl = interm_lvl,
  device = default_device()
)

# -----------------------------
# STEP 4: EXCESS PROBABILITY
# -----------------------------

# empirical 0.99 quantile of test observations
emp_quantile <- as.numeric(quantile(y_test, probs = 0.99, na.rm = TRUE))

ppred_eqrn <- EQRN_excess_probability(
  val = emp_quantile,
  fit_eqrn = fit_EQRN,
  X = X_test,
  intermediate_quantiles = intermediateq_test,
  interm_lvl = interm_lvl,
  proba_type = "excess"
)

# Check
i <- which.max(y_test)
cat("Peak gust:", y_test[i])
cat("Predicted exceedance prob:", ppred_eqrn[i], "\n")
cat("Predicted 0.99 quantile:", qpred_eqrn[i], "\n")

# ------------------------------------------------------------
# Compute baseline using quantile regression forest
# ------------------------------------------------------------
# Intermediate quantile model prediction of extreme quantile
qpred_extreme_grf  <- predict(fit_grf, newdata = X_test, quantiles = levels_predict)$predictions

