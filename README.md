# Extreme Quantile Regression for Wind Gust Forecasts (CIENS + EQRN)

This repository contains the code for a seminar project applying **Extreme Quantile Regression Networks (EQRN)** to the prediction of **extreme wind gusts** using the **CIENS ensemble forecast dataset**.

The goal is to model conditional extreme quantiles (in particular the **0.99 quantile**) of wind gust speeds and to evaluate calibration and tail behavior across different German stations.

---

## Project Overview

The workflow follows the two-stage EQRN approach:

1. **Intermediate Quantile Model (Model 1)**  
   Estimate a moderate quantile level (e.g. τ₀ = 0.8) using:
   - Generalized Random Forest (IID setup)
   - Quantile Regression Network (sequential setup)

2. **Tail Model (Model 2)**  
   Fit a Generalized Pareto tail model via EQRN to obtain extreme quantiles for τ > τ₀.

The project compares:

- IID EQRN vs GRF baseline  
- Sequential (RNN) EQRN vs QRN baseline  
- Empirical test quantile and ensemble maximum as reference thresholds  

---

## Requirements

Packages are managed using **renv**.

To restore the required environment, run the following in the R console:

```r
install.packages("renv")
renv::restore()
```

## Reproducing Results

To reproduce the results, run the following scripts in order:

- `data_download.R`  
- `scripts/data_preprocessing.R`  
- Model training:  
  - IID setup: `scripts/GRF_EQRN_train.R`  
  - Sequential setup: `scripts/QRN_EQRN_train.R`  
- Model evaluation:  
  - IID setup: `scripts/evaluation.R`  
  - Sequential setup: `scripts/evaluation_seq.R`  
- Summary table generation:  
  - `scripts/coverage_table.R`  

---

## References

Lerch, S., Schulz, B., Hess, R., Moeller, A., Primo, C., Trepte, S., & Theis, S. (2026).  
*Operational convection-permitting COSMO/ICON ensemble predictions at observation sites (CIENS).*  
Geoscience Data Journal, 13(1), e70051.

Pasche, O., & Engelke, S. (2024).  
*Neural networks for extreme quantile regression with an application to forecasting of flood risk.*  
The Annals of Applied Statistics, 18(4), 2818–2839.

