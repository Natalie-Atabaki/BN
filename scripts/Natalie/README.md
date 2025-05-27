# Bayesian Network Analysis of Clinical and Proteomic Data

This repository contains two main R scripts for building and analysing Bayesian Networks (BNs) using clinical and proteomic data.

---

## 📁 Files

### 1. `Building_Networks.R`
This script performs the following:
- Loads and merges clinical and protein datasets
- Fits Bayesian Networks using the `bnlearn` package
- Performs k-fold cross-validation with parallel computing
- Computes network structure strength and arc direction probabilities
- Filters networks based on edge strength (e.g., ≥0.8) and directionality
- Saves networks and arc strengths to `.RDS` and `.xlsx` files

**Key Outputs:**
- `network.rds`: averaged network and its strengths
- `ArcStrengthDirection.xlsx`: arc strengths at various thresholds

---

### 2. `Posterior_Prob.R`
This script focuses on:
- Simulating or importing a discretised clinical-protein dataset
- Fitting a BN to the discrete data using Bayesian parameter estimation
- Estimating posterior probabilities of a trait (e.g., `Trait1`- `LiverFat' in the paper ) given levels (HIGH/LOW) of selected proteins
- Visualising posterior probabilities as a polar bar chart

**Key Features:**
- Uses `gRain` for inference
- Divides proteins into `High` and `Low` based on hypothesised direction of association
- Highlights probabilistic relationships between traits and individual proteins

---

## 📦 Dependencies

Both scripts require the following R packages:

```r
install.packages(c(
  "bnlearn", "gRain", "ggplot2", "dplyr", 
  "openxlsx", "parallel", "Rgraphviz"
))

