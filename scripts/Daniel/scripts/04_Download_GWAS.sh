#!/bin/bash

DIR=/Users/da1078co/Documents/Lund/PhD/Projects/BN
SUMDIR=${DIR}/data/GWAS_SumStat

TRAIT=BasalISR
# Basal ISR from IVGTT GWAS from Wood et al 2017
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://gwas.mrcieu.ac.uk/files/ebi-a-GCST004488/ebi-a-GCST004488.vcf.gz

TRAIT=Glucose
# Fasting insulin from MAGIC consortium Chen et al 2021
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002232/GCST90002232_buildGRCh37.tsv.gz

TRAIT=HDL
# HDL from GLGC consortium Graham et al 2021 
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/trans_ancestry/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_HDL_INV_ALL_with_N_1.gz

TRAIT=HOMA_IR
# HOMA-IR from MAGIC consortium Dupuis et al 2011
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://gwas.mrcieu.ac.uk/files/ieu-b-118/ieu-b-118.vcf.gz

TRAIT=HbA1c
# HbA1c from MAGIC consortium Chen et al 2021
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002244/GCST90002244_buildGRCh37.tsv.gz

TRAIT=Insulin
# Fasting insulin from MAGIC consortium Chen et al 2021
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002238/GCST90002238_buildGRCh37.tsv.gz

TRAIT=LiverFat
# Liver fat from UKBB GWAS from Liu et al 2021
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90016001-GCST90017000/GCST90016673/GCST90016673_buildGRCh37.tsv.gz

TRAIT=PancFat
# Pancreatic fat from UKBB GWAS from Liu et al 2021
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90016001-GCST90017000/GCST90016675/GCST90016675_buildGRCh37.tsv.gz

TRAIT=SAT
# Subcutaneous adipose tissue from UKBB GWAS from Liu et al 2021
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90016001-GCST90017000/GCST90016672/GCST90016672_buildGRCh37.tsv.gz

TRAIT=TG
# TG from GLGC consortium Graham et al 2021
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/trans_ancestry/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_logTG_INV_ALL_with_N_1.gz

TRAIT=TwoGlucose
# 2 hour glucose from MAGIC consortium Chen et al 2021
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002227/GCST90002227_buildGRCh37.tsv.gz

TRAIT=TwoInsulin
# 2 hour insulin fold change from MAGIC consortium Chen et al 2021
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://magicinvestigators.org/downloads/MAGIC_postchallengeIR_IFC_noBMI_ALL.tsv.gz

TRAIT=VAT
# Visceral adipose tissue from UKBB GWAS from Liu et al 2021
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90016001-GCST90017000/GCST90016671/GCST90016671_buildGRCh37.tsv.gz

TRAIT=BETACELL
# Beta cell function GWAS from Madsen et al 2024
# From this we will extract xinsdG30, equivalent to GlucoseSens
mkdir $SUMDIR/$TRAIT
curl -L \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://api.kpndataregistry.org/api/d/ChTEgN










