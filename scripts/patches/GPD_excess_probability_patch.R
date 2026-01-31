# Patch so the excess probability at exceedances is not a character

GPD_excess_probability <- function(val, sigma, xi, interm_threshold, threshold_p,
                                   body_proba = "default",
                                   proba_type = c("excess","cdf")) {
  
  proba_type <- match.arg(proba_type)
  
  n <- max(length(val), length(sigma), length(xi), length(interm_threshold))
  for(v in list(val, sigma, xi, interm_threshold)){
    if(length(v) != 1 && length(v) != n)
      stop("val, sigma, xi, interm_threshold must be equal or unit length.")
  }
  val <- rep_len(c(val), n)
  sigma <- rep_len(c(sigma), n)
  xi <- rep_len(c(xi), n)
  interm_threshold <- rep_len(c(interm_threshold), n)
  
  Prob <- rep(NA_real_, n)
  
  # exact threshold point
  Prob[val == interm_threshold] <- threshold_p
  
  # tail part: val > threshold
  inds_ex <- val > interm_threshold
  Z <- val - interm_threshold
  gpd_cdf <- mapply(evd::pgpd, q = Z[inds_ex], loc = 0,
                    scale = sigma[inds_ex], shape = xi[inds_ex], SIMPLIFY = TRUE)
  Prob[inds_ex] <- threshold_p + gpd_cdf * (1 - threshold_p)
  
  if (proba_type == "excess") {
    Prob[val >= interm_threshold] <- 1 - Prob[val >= interm_threshold]
  }
  
  # body part: val < threshold
  if (any(val < interm_threshold)) {
    if (is.character(body_proba) && body_proba == "default") {
      # numeric default instead of strings like ">0.2"
      Prob[val < interm_threshold] <- (1 - threshold_p)  # for excess
      if (proba_type == "cdf") Prob[val < interm_threshold] <- threshold_p
    } else if (is.numeric(body_proba)) {
      Prob[val < interm_threshold] <- body_proba
    } else {
      Prob[val < interm_threshold] <- NA_real_
    }
  }
  
  Prob
}

EQRN_excess_probability <- function(val, fit_eqrn, X, intermediate_quantiles,
                                    interm_lvl = fit_eqrn$interm_lvl,
                                    body_proba = "default",
                                    proba_type = c("excess","cdf"),
                                    device = default_device()) {
  
  proba_type <- match.arg(proba_type)
  
  # --- ensure body_proba is numeric (avoid strings like ">0.2") -------------
  if (is.character(body_proba) && body_proba == "default") {
    body_proba <- if (proba_type == "excess") (1 - interm_lvl) else interm_lvl
  }
  # If user passes a string like ">0.2", strip it and convert to numeric
  if (is.character(body_proba)) {
    body_proba <- suppressWarnings(as.numeric(gsub("^\\s*[<>]=?\\s*", "", body_proba)))
  }
  if (!is.numeric(body_proba) || length(body_proba) != 1L || !is.finite(body_proba)) {
    stop("body_proba must be numeric (or 'default').")
  }
  
  GPD_params_pred <- EQRN_predict_params(
    fit_eqrn, X, intermediate_quantiles,
    return_parametrization = "classical",
    interm_lvl = interm_lvl,
    device = device
  )
  sigmas <- as.numeric(GPD_params_pred$scales)
  xis    <- as.numeric(GPD_params_pred$shapes)
  
  Probs <- GPD_excess_probability(
    val,
    sigma = sigmas,
    xi = xis,
    interm_threshold = intermediate_quantiles,
    threshold_p = interm_lvl,
    body_proba = body_proba,      # now guaranteed numeric
    proba_type = proba_type
  )
  
  # guarantee numeric vector (in case upstream coerced)
  Probs <- suppressWarnings(as.numeric(Probs))
  return(Probs)
}

EQRN_excess_probability_seq <- function(val, fit_eqrn, X, Y, intermediate_quantiles,
                                        interm_lvl = fit_eqrn$interm_lvl,
                                        crop_predictions = FALSE,
                                        body_proba = "default",
                                        proba_type = c("excess","cdf"),
                                        seq_len = fit_eqrn$seq_len,
                                        device = default_device(),
                                        eps_sigma = 1e-6,
                                        na_to_zero_when_out_of_support = TRUE) {
  
  proba_type <- match.arg(proba_type)
  
  GPD_params_pred <- EQRN_predict_params_seq(
    fit_eqrn, X, Y, intermediate_quantiles,
    return_parametrization = "classical",
    interm_lvl = interm_lvl,
    seq_len = seq_len,
    device = device
  )
  
  sigmas <- as.numeric(GPD_params_pred$scales)
  xis    <- as.numeric(GPD_params_pred$shapes)
  
  # --- align thresholds to the prediction vector length --------------------
  # GPD_excess_probability expects one threshold per predicted probability (after cropping)
  u <- as.numeric(intermediate_quantiles)
  u_tail <- u[(seq_len + 1):length(u)]  # same as your original
  
  # --- basic parameter sanitizing ------------------------------------------
  # enforce sigma > 0
  sigmas[!is.finite(sigmas)] <- NA_real_
  sigmas <- pmax(sigmas, eps_sigma)
  
  # if xi is non-finite, set to 0 (GPD limit case); alternatively set NA
  xis[!is.finite(xis)] <- 0
  
  Probs <- GPD_excess_probability(
    val,
    sigma = sigmas,
    xi = xis,
    interm_threshold = u_tail,
    threshold_p = interm_lvl,
    body_proba = body_proba,
    proba_type = proba_type
  )
  
  # --- fix domain failures for bounded tails (xi < 0) -----------------------
  # If xi < 0, support is z < -sigma/xi, i.e. val < u - sigma/xi.
  # When val is beyond that endpoint, the exceedance probability should be 0.
  if (na_to_zero_when_out_of_support) {
    z_val <- val - u_tail
    out_of_support <- (xis < 0) & is.finite(z_val) & (z_val >= (-sigmas / xis))
    
    Probs[out_of_support] <- 0
  }
  
  # --- crop back to original length if needed ------------------------------
  if (!crop_predictions) {
    Probs <- c(rep(NA_real_, seq_len), Probs)
  }
  
  Probs
}

