# Extreme Quantile Regression for Wind Gust Forecasts (CIENS + EQRN)

This repository contains the code for a seminar project applying **Extreme Quantile Regression Networks (EQRN)** to the prediction of **extreme wind gusts** using the **CIENS ensemble forecast dataset**.

The goal is to model conditional extreme quantiles (in particular the **0.99 quantile**) of wind gust speeds and to evaluate calibration and tail behavior across different German stations.

---

The workflow follows the two-stage EQRN approach:

1. **Intermediate Quantile Model (Model 1)**  
   Estimate a moderate quantile level (e.g. τ₀ = 0.8) using:
   - Generalized Random Forest (IID setup)
   - Quantile Regression Network (Sequential setup)

2. **Tail Model (Model 2)**  
   Fit a Generalized Pareto tail model via EQRN to obtain extreme quantiles for τ > τ₀.

The project compares:

- IID EQRN vs GRF baseline  
- Sequential (RNN) EQRN vs QRN baseline  
- Empirical test quantile and ensemble maximum as reference thresholds

