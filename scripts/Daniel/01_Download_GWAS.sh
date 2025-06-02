#!/bin/bash

DIR=~/
SUMDIR=${DIR}/data/GWAS_SumStat

TRAIT=AbdSAT
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90016001-GCST90017000/GCST90016672/GCST90016672_buildGRCh37.tsv.gz

TRAIT=AbdVAT
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90016001-GCST90017000/GCST90016671/GCST90016671_buildGRCh37.tsv.gz

TRAIT=BasalISR
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://gwas.mrcieu.ac.uk/files/ebi-a-GCST004488/ebi-a-GCST004488.vcf.gz

TRAIT=BETACELL
mkdir $SUMDIR/$TRAIT
curl -L \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://api.kpndataregistry.org/api/d/ChTEgN

TRAIT=FG
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002232/GCST90002232_buildGRCh37.tsv.gz

TRAIT=FI
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002238/GCST90002238_buildGRCh37.tsv.gz

TRAIT=Glu2rh
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002227/GCST90002227_buildGRCh37.tsv.gz

TRAIT=HbA1c
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002244/GCST90002244_buildGRCh37.tsv.gz

TRAIT=HDL
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/trans_ancestry/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_HDL_INV_ALL_with_N_1.gz

TRAIT=HOMA_IR
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://gwas.mrcieu.ac.uk/files/ieu-b-118/ieu-b-118.vcf.gz

TRAIT=IFC
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://magicinvestigators.org/downloads/MAGIC_postchallengeIR_IFC_noBMI_ALL.tsv.gz


TRAIT=ISI
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://magicinvestigators.org/downloads/MAGIC_postchallengeIR_ISI_noBMI_ALL.tsv.gz


TRAIT=LiverFat
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90016001-GCST90017000/GCST90016673/GCST90016673_buildGRCh37.tsv.gz


TRAIT=PancreasFat
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90016001-GCST90017000/GCST90016675/GCST90016675_buildGRCh37.tsv.gz


TRAIT=TG
mkdir $SUMDIR/$TRAIT
curl \
--remote-name \
--output-dir $SUMDIR/$TRAIT \
https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/trans_ancestry/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_logTG_INV_ALL_with_N_1.gz




