#!/bin/bash

# Download reference panel
# This is derived from 1KG removing indels and only MAF > 1%
curl http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz --output /Users/da1078co/RefData/1kg.v3.tgz

# Uncompress
tar -xvzf /Users/da1078co/RefData/1kg.v3.tgz

# Calculate allele frequencies
/Users/da1078co/miniforge3/pkgs/plink2-2.00a5.12-hac4f329_0/bin/plink2 \
--bfile /Users/da1078co/RefData/1KG/EUR \
--freq cols=chrom,pos,ref,alt,altfreq \
--out /Users/da1078co/RefData/1KG/EUR

# Select nonpalindromic SNPs
awk \
'\
 BEGIN{OFS="\t"}\
 !(\
   (\
    ($3=="A" && $4=="T") || \
    ($3=="T" && $4=="A") || \
    ($3=="C" && $4=="G") || \
    ($3=="G" && $4=="C") \
   ) &&\
   ($6 >= 0.4 || $6 <= 0.6)\
 )\
' \
/Users/da1078co/RefData/1KG/EUR.afreq > \
/Users/da1078co/RefData/1KG/EUR_NOPAL.tsv

# Download GENCODE mappings

## Basic
curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh37_mapping/gencode.v47lift37.basic.annotation.gtf.gz --output /Users/da1078co/Documents/Data/GENCODE/gencode.v47lift37.basic.annotation.gtf.gz

## Extended
curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh37_mapping/gencode.v47lift37.basic.annotation.gff3.gz --output /Users/da1078co/Documents/Data/GENCODE/gencode.v47lift37.basic.annotation.gff3.gz