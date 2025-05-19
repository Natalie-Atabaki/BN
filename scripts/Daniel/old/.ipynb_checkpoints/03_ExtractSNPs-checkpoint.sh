#!/bin/bash

DIR=/Users/da1078co/Documents/PhD/Projects/BN_Naeimeh
SUMDIR=${DIR}/data/GWAS_SumStat
RESDIR=${DIR}/data/GWAS_Hits
REFPANEL=/Users/da1078co/RefData/1KG/EUR_NOPAL.tsv
HEADER="CHROM\tPOS\tEA\tNEA\tEAF\tBETA\tSE\tPVAL\tRSID\tREFA1\tREFA2\tREFA2FREQ"

TRAIT=AbdSAT
SUMSTAT=GCST90016672_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hits.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk '$2 <= 1e-6' |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$3"\t"$4} \
(b in a){print $3,$4,$5,$6,$7,$8,$9,$2,a[b]}' $REFPANEL - >> $RES

TRAIT=AbdVAT
SUMSTAT=GCST90016671_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hits.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk '$2 <= 1e-6' |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$3"\t"$4} \
(b in a){print $3,$4,$5,$6,$7,$8,$9,$2,a[b]}' $REFPANEL - >>\
$RES

TRAIT=BasalISR
SUMSTAT=ebi-a-GCST004488.vcf.gz
RES=${RESDIR}/${TRAIT}_hits.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
grep -v "#" |\
awk 'BEGIN{OFS="\t"} {print $1,$2,$4,$5,$3,$10}' |\
sed 's/:/\t/g' |\
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$9,$6,$7,10**(-$8)}' |\
awk '$8 <= 1e-6' |\
awk 'BEGIN{OFS="\t"}\
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$1"\t"$2}\
(b in a){print $0,a[b]}' $REFPANEL - >>\
$RES

TRAIT=FG
SUMSTAT=GCST90002232_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hits.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk '$8 <= 1e-6' |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$1"\t"$2} \
(b in a){print $1,$2,$3,$4,$5,$6,$7,$8,a[b]}' $REFPANEL - >>\
$RES

TRAIT=FI
SUMSTAT=GCST90002238_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hits.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk '$8 <= 1e-6' |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$1"\t"$2} \
(b in a){print $1,$2,$3,$4,$5,$6,$7,$8,a[b]}' $REFPANEL - >>\
$RES

TRAIT=GLCGN
SUMSTAT=Glucagon_4891_50.txt.gz
RES=${RESDIR}/${TRAIT}_hits.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk '$11 <= 1e-6' |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$18"\t"$19} \
(b in a){print $18,$19,toupper($3),toupper($4),$5,$9,$10,$11,a[b]}' $REFPANEL - >>\
$RES

TRAIT=GLP1
SUMSTAT=GLP1R_13085_18.txt.gz
RES=${RESDIR}/${TRAIT}_hits.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk '$11 <= 1e-6' |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$18"\t"$19} \
(b in a){print $18,$19,toupper($3),toupper($4),$5,$9,$10,$11,a[b]}' $REFPANEL - >>\
$RES

TRAIT=Glu2hr
SUMSTAT=GCST90002227_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hits.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk '$8 <= 1e-6' |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$1"\t"$2} \
(b in a){print $1,$2,$3,$4,$5,$6,$7,$8,a[b]}' $REFPANEL - >>\
$RES

TRAIT=HbA1c
SUMSTAT=GCST90002244_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hits.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk '$8 <= 1e-6' |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$1"\t"$2} \
(b in a){print $1,$2,$3,$4,$5,$6,$7,$8,a[b]}' $REFPANEL - >>\
$RES

TRAIT=HDL
SUMSTAT=with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_HDL_INV_ALL_with_N_1.gz
RES=${RESDIR}/${TRAIT}_hits.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk '$16 <= 1e-6' |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$2"\t"$3} \
(b in a){print $2,$3,$4,$5,$8,$14,$15,$16,a[b]}' $REFPANEL - >>\
$RES

TRAIT=HOMA_IR
SUMSTAT=ieu-b-118.vcf.gz
RES=${RESDIR}/${TRAIT}_hits.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
grep -v "#" |\
awk 'BEGIN{OFS="\t"} {print $1,$2,$4,$5,$3,$10}' |\
sed 's/:/\t/g' |\
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,".",$6,$7,10**(-$8)}' |\
awk '$8 <= 1e-6' |\
awk 'BEGIN{OFS="\t"}\
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$1"\t"$2}\
(b in a){print $0,a[b]}' $REFPANEL - >>\
$RES

TRAIT=IFC
SUMSTAT=MAGIC_postchallengeIR_IFC_noBMI_ALL.tsv.gz
RES=${RESDIR}/${TRAIT}_hits.tsv
echo -e $HEADER > $RES
awk '$8<=1e-6' ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"}\
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$1"\t"$2}\
(b in a){print $1,$2,$3,$4,$7,$5,$6,$8,a[b]}' $REFPANEL - >>\
$RES

TRAIT=ISI
SUMSTAT=MAGIC_postchallengeIR_ISI_noBMI_ALL.tsv.gz
RES=${RESDIR}/${TRAIT}_hits.tsv
echo -e $HEADER > $RES
awk '$8<=1e-6' ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"}\
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$1"\t"$2}\
(b in a){print $1,$2,$3,$4,$7,$5,$6,$8,a[b]}' $REFPANEL - >>\
$RES

TRAIT=LiverFat
SUMSTAT=GCST90016673_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hits.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk '$2 <= 1e-6' |\
awk 'BEGIN{OFS="\t"}\
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$3"\t"$4}\
(b in a){print $3,$4,$5,$6,$7,$8,$9,$2,a[b]}' $REFPANEL - >>\
$RES

TRAIT=PancreasFat
SUMSTAT=GCST90016675_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hits.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk '$2 <= 1e-6' |\
awk 'BEGIN{OFS="\t"}\
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$3"\t"$4}\
(b in a){print $3,$4,$5,$6,$7,$8,$9,$2,a[b]}' $REFPANEL - >>\
$RES

TRAIT=TG
SUMSTAT=with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_logTG_INV_ALL_with_N_1.gz
RES=${RESDIR}/${TRAIT}_hits.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk '$16 <= 1e-6' |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$2"\t"$3} \
(b in a){print $2,$3,$4,$5,$8,$14,$15,$16,a[b]}' $REFPANEL - >>\
$RES