# Bayesian Network Analysis of Clinical and Proteomic Data

This repository contains two main R scripts for building and analysing Bayesian Networks (BNs) using simulated clinical and proteomic data.

---

## 📁 Files

### 1. `Building_Networks.R`
This script performs the following:
- Simulates and merges clinical and proteomic datasets
- Fits Bayesian Networks using the `bnlearn` package
- Conducts k-fold cross-validation with parallel computing
- Estimates arc strengths and direction probabilities from ensemble models
- Filters networks based on confidence thresholds (strength ≥ 0.8, direction ≥ 0.5 or 0.8)
- Compares learned network BIC scores to randomly generated networks
- Generates visual output for BIC score distribution

---

### 2. `Posterior_Prob.R`
This script focuses on:
- Simulating a discretised clinical-protein dataset (selected features from the previous step)
- Fitting a discrete BN to the data using Bayesian parameter estimation
- Computing posterior probabilities of a selected trait (e.g., `Trait1`)
  given HIGH or LOW expression levels of selected proteins
- Visualising the inferred probabilities with a polar bar chart

---

## 📦 Dependencies

Both scripts rely on the following R packages:

```r
install.packages(c(
  "bnlearn", "gRain", "ggplot2", "dplyr", 
  "openxlsx", "parallel"
))
