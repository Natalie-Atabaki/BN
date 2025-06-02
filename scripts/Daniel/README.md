# MR analysis

This part of the repository contains the scripts used to test the links suggested by the Bayesian Network using summary data.

## Files

|File                        |Description                                                   |
|----------------------------|--------------------------------------------------------------|
|`00_RefMaps.sh`             |Download of reference data (e.g. reference panel              |
|`01_Download_GWAS.sh`       |Download of GWAS summary statistics                           |
|`02_Download_pQTLPanel.py`  |Download of pQTL meta information                             |
|`03_Process_pQTLPanel.ipynb`|Process pQTL meta information                                 |
|`04_Download_pQTL.py`       |Download of pQTL summary data                                 |
|`05_Harmonize.sh`           |Harmonizing GWAS summary statistics                           |
|`06_UTILS.R`                |Functions used to run MR                                      |
|`07_PrepareEdges.ipynb`     |Preparing edges from the Bayesian Network to be tested        |
|`08_EdgeTest.R`             |Running MR analyses on selected edges of the Bayesian Network |
|`09_EdgeResults.ipynb`      |Processing results                                            |

## Dependencies

- `GNU bash, version 5.1.8(1)-release (x86_64-redhat-linux-gnu)`

- 
