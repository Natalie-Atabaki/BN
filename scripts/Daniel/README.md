# MR analysis

This part of the repository contains the scripts used to test the links suggested by the Bayesian Network using summary data.

## Files

|File                         |Description                                                   |
|-----------------------------|--------------------------------------------------------------|
|`00_Functions.R`             |Functions used in GWAS harmonization and MR                   |
|`01_PrepareEdges.ipynb`      |Preparing edges from the Bayesian Network to be tested        |
|`02_RefMaps.sh`              |Download of reference data (e.g. reference panel)             |
|`03_Download_GWAS.sh`        |Download of GWAS summary statistics                           |
|`04a_Download_pQTLPanel.py`  |Download of pQTL meta information                             |
|`04b_Process_pQTLPanel.ipynb`|Process pQTL meta information                                 |
|`04c_Download_pQTL.py`       |Download of pQTL summary data                                 |
|`05_Harmonize.sh`            |Harmonizing GWAS summary statistics                           |
|`06_EdgeTest.R`              |Running MR analyses on selected edges of the Bayesian Network |
|`07_EdgeResults.ipynb`       |Processing results                                            |

## Dependencies

- `GNU bash, version 3.2.57(1)-release (arm64-apple-darwin24)`

- `PLINK v2.00a5.12 64-bit (25 Jun 2024)`

- `R version 4.4.3 (2025-02-28) -- "Trophy Case"`
  - `here==1.0.1`
  - `readr==2.1.5`
  - `dplyr==1.1.4`
  - `tidyr==1.3.1`
  - `purrr==1.0.2`
  - `survey==4.4.2`
  - `ggplot2==3.5.1`
  - `stringr==1.5.1`

- `Python 3.12.5`:
  - `synapseclient==4.6.0`
  - `pandas==2.2.3`
  - `tarfile=0.9.0`
