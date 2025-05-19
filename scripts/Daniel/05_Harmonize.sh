#!/bin/bash

DIR=/Users/da1078co/Documents/Lund/PhD/Projects/BN_Naeimeh
SUMDIR=${DIR}/data/GWAS_SumStat
RESDIR=${DIR}/data/GWAS_Harmonized
REFPANEL=/Users/da1078co/Documents/Data/1KG/EUR_NOPAL.tsv
HEADER="CHROM\tPOS\tEA\tNEA\tEAF\tBETA\tSE\tPVAL\tRSID"

HarmonFx() {
  awk \
  'BEGIN{OFS="\t"} \
  {
    CHROM = $1
    POS = $2
    EA = $3
    NEA = $4
    EAF = $5
    BETA = $6
    SE = $7
    PVAL = $8
    RSID = $9
    REFA1 = $10
    REFA2 = $11
    REFA2FREQ = $12
    
    if (BETA > 0) {
      EA_I = EA
      NEA_I = NEA
      EAF_I = EAF
    } else {
      EA_I = NEA
      NEA_I = EA
      EAF_I = 1 - EAF
    }

    BETA_I = (BETA < 0) ? -BETA : BETA

    PALIND1 = (REFA1 == "A" && REFA2 == "T") || (REFA1 == "T" && REFA2 == "A")
    PALIND2 = (REFA1 == "C" && REFA2 == "G") || (REFA1 == "G" && REFA2 == "C")
    PALIND = PALIND1 || PALIND2

    REFA1C = (REFA1 == "A") ? "T": (REFA1 == "T") ? "A" :(REFA1 == "C") ? "G" : "C"
    REFA2C = (REFA2 == "A") ? "T": (REFA2 == "T") ? "A" :(REFA2 == "C") ? "G" : "C"

    REFA2FREQF = 1 - REFA2FREQ
    
    EAFDIFF1 = REFA2FREQ - EAF_I
    EAFDIFF1 = sqrt(EAFDIFF1*EAFDIFF1)
    
    EAFDIFF2 = REFA2FREQF - EAF_I
    EAFDIFF2 = sqrt(EAFDIFF2*EAFDIFF2)

    if ( EA_I == REFA2 && NEA_I == REFA1 ) {
      if ( !(EAF ~ /^0\\.[0-9]+$/) ) {
        if ( !PALIND ) {
          REFA2FREQH = REFA2FREQ
        } else { next }
      } else if ( EAFDIFF1 < 0.2 ) {
        REFA2FREQH = REFA2FREQ
      } else if ( PALIND && EAFDIFF2 < 0.2 ) {
        REFA2FREQH = REFA2FREQF
      } else { next } 
    } else if ( EA_I == REFA1 && NEA_I == REFA2 ) {
      if ( !(EAF ~ /^0\\.[0-9]+$/) ) {
        if ( !PALIND ) {
          REFA2FREQH = REFA2FREQF
        } else { next }
      } else if ( EAFDIFF2 < 0.2 ) {
        REFA2FREQH = REFA2FREQF
      } else if ( PALIND && EAFDIFF1 < 0.2 ){
        REFA2FREQH = REFA2FREQ
      } else { next } 
    } else if ( !PALIND ){
      if ( EA_I == REFA2C && NEA_I == REFA1C ){
        REFA2FREQH = REFA2FREQ
      } else if ( EA_I == REFA1C && NEA_I == REFA2C ){
        REFA2FREQH = REFA2FREQF
      } else { next } 
    } else { next } 
    
    EAFX = ( !(EAF ~ /^0\\.[0-9]+$/) ) ? REFA2FREQH : EAF_I

    EAFDIFF = REFA2FREQH - EAFX
    EAFDIFF = sqrt(EAFDIFF*EAFDIFF)

    if( EAFDIFF < 0.2 ){
      print CHROM, POS, EA_I, NEA_I, EAFX, BETA_I, SE, PVAL, RSID
    }
    
  }' \
  $1
}

TRAIT=AbdSAT
SUMSTAT=GCST90016672_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$3"\t"$4} \
(b in a){print $3,$4,$5,$6,$7,$8,$9,$2,a[b]}' $REFPANEL - |\
HarmonFx |\
sort -k1,1V -k2,2n | uniq >> $RES

TRAIT=AbdVAT
SUMSTAT=GCST90016671_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$3"\t"$4} \
(b in a){print $3,$4,$5,$6,$7,$8,$9,$2,a[b]}' $REFPANEL - |\
HarmonFx |\
sort -k1,1V -k2,2n | uniq >> $RES

TRAIT=BasalISR
SUMSTAT=ebi-a-GCST004488.vcf.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
grep -v "#" |\
awk 'BEGIN{OFS="\t"} {print $1,$2,$5,$4,$10}' |\
sed 's/:/\t/g' |\
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$8,$5,$6,10**(-$7)}' |\
awk 'BEGIN{OFS="\t"}\
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$1"\t"$2}\
(b in a){print $0,a[b]}' $REFPANEL - |\
HarmonFx |\
sort -k1,1V -k2,2n | uniq >> $RES

TRAIT=CIR
SUMSTAT=ChTEgN
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
unzip -c $SUMDIR/BETACELL/${SUMSTAT} meta1cir.filthetmafn.rsid.selectedcolumns |\
tail -n+4 |\
sed 's/:/\t/g' |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$1"\t"$2} \
(b in a){print $1,$2,toupper($3),toupper($4),$5,$9,$10,$11,a[b]}' $REFPANEL - |\
HarmonFx |\
sort -k1,1V -k2,2n | uniq >> $RES

TRAIT=FG
SUMSTAT=GCST90002232_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$1"\t"$2} \
(b in a){print $1,$2,$3,$4,$5,$6,$7,$8,a[b]}' $REFPANEL - |\
HarmonFx |\
sort -k1,1V -k2,2n | uniq >> $RES

TRAIT=FI
SUMSTAT=GCST90002238_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$1"\t"$2} \
(b in a){print $1,$2,$3,$4,$5,$6,$7,$8,a[b]}' $REFPANEL - |\
HarmonFx |\
sort -k1,1V -k2,2n | uniq >> $RES

TRAIT=Glu2hr
SUMSTAT=GCST90002227_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$1"\t"$2} \
(b in a){print $1,$2,$3,$4,$5,$6,$7,$8,a[b]}' $REFPANEL - |\
HarmonFx |\
sort -k1,1V -k2,2n | uniq >> $RES

TRAIT=HbA1c
SUMSTAT=GCST90002244_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$1"\t"$2} \
(b in a){print $1,$2,$3,$4,$5,$6,$7,$8,a[b]}' $REFPANEL - |\
HarmonFx |\
sort -k1,1V -k2,2n | uniq >> $RES

TRAIT=HDL
SUMSTAT=with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_HDL_INV_ALL_with_N_1.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$2"\t"$3} \
(b in a){print $2,$3,$5,$4,$8,$14,$15,$16,a[b]}' $REFPANEL - |\
HarmonFx |\
sort -k1,1V -k2,2n | uniq >> $RES

TRAIT=HOMA_IR
SUMSTAT=ieu-b-118.vcf.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
grep -v "#" |\
awk 'BEGIN{OFS="\t"} {print $1,$2,$5,$4,$10}' |\
sed 's/:/\t/g' |\
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,".",$5,$6,10**(-$7)}' |\
awk 'BEGIN{OFS="\t"}\
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$1"\t"$2}\
(b in a){print $0,a[b]}' $REFPANEL - |\
HarmonFx |\
sort -k1,1V -k2,2n | uniq >> $RES

TRAIT=IFC
SUMSTAT=MAGIC_postchallengeIR_IFC_noBMI_ALL.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
awk 'BEGIN{OFS="\t"}\
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$1"\t"$2}\
(b in a){print $1,$2,$3,$4,$7,$5,$6,$8,a[b]}' $REFPANEL ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
HarmonFx |\
sort -k1,1V -k2,2n | uniq >> $RES

TRAIT=ISI
SUMSTAT=MAGIC_postchallengeIR_ISI_noBMI_ALL.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
awk 'BEGIN{OFS="\t"}\
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$1"\t"$2}\
(b in a){print $1,$2,$3,$4,$7,$5,$6,$8,a[b]}' $REFPANEL ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
HarmonFx |\
sort -k1,1V -k2,2n | uniq >> $RES

TRAIT=LiverFat
SUMSTAT=GCST90016673_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"}\
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$3"\t"$4}\
(b in a){print $3,$4,$5,$6,$7,$8,$9,$2,a[b]}' $REFPANEL - |\
HarmonFx |\
sort -k1,1V -k2,2n | uniq >> $RES

TRAIT=PancreasFat
SUMSTAT=GCST90016675_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"}\
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$3"\t"$4}\
(b in a){print $3,$4,$5,$6,$7,$8,$9,$2,a[b]}' $REFPANEL - |\
HarmonFx |\
sort -k1,1V -k2,2n | uniq >> $RES

TRAIT=TG
SUMSTAT=with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_logTG_INV_ALL_with_N_1.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$2"\t"$3} \
(b in a){print $2,$3,$5,$4,$8,$14,$15,$16,a[b]}' $REFPANEL - |\
HarmonFx |\
sort -k1,1V -k2,2n | uniq >> $RES

TRAIT=xinsdG30
SUMSTAT=ChTEgN
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
unzip -c $SUMDIR/BETACELL/${SUMSTAT} meta1xinsdG30.filthetmafn.rsid.selectedcolumns |\
tail -n+4 |\
sed 's/:/\t/g' |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{b=$1"\t"$2} \
(b in a){print $1,$2,toupper($3),toupper($4),$5,$9,$10,$11,a[b]}' $REFPANEL - |\
HarmonFx |\
sort -k1,1V -k2,2n | uniq >> $RES

PQTLPANEL=${DIR}/data/protquery_olinkUKB.tsv

tail -n+2 $PQTLPANEL |\
cut -f1,18 |\
sed 's/.tar//g' |\
while read -r GENE GFOLDER
do
   RES=${RESDIR}/${GENE}_hrmnzd.tsv
   echo -e $HEADER > $RES
   for CHROMFILE in ${SUMDIR}/${GENE}/${GFOLDER}/*
   do
     gunzip -c ${CHROMFILE} |\
     tail -n+2 |\
     sed -E 's/[A-Z0-9]+:([0-9]+):[a-zA-Z0-9:]+ /\1 /g' 
   done |\
   awk \
   '(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
   {b=$1"\t"$3} \
   (b in a){print $1,$3,$5,$4,$6,$10,$11,10**(-$13),a[b]}' \
   $REFPANEL - |\
   HarmonFx |\
   sort -k1,1V -k2,2n | uniq >> $RES
done