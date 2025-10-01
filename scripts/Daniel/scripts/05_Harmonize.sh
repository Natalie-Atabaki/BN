#!/bin/bash

# Directories and files
DIR=/Users/da1078co/Documents/Lund/PhD/Projects/BN
SUMDIR=${DIR}/data/GWAS_SumStat
RESDIR=${DIR}/data/GWAS_Harmonized
REFPANEL=/Users/da1078co/Documents/Data/1KG/EUR_NOPAL.tsv
HEADER="CHROM\tPOS\tEA\tNEA\tEAF\tBETA\tSE\tPVAL\tRSID"

# Harmonization function
HarmonFx() {
  awk \
  'BEGIN {
    OFS="\t"
  }
  {
    # Expected input columns
    CHROM = $1; POS = $2
    EA = $3; NEA = $4; EAF = $5
    BETA = $6; SE = $7; PVAL = $8
    RSID = $9
    REFA1 = $10; REFA2 = $11; REFA2FREQ = $12
    
    # Aligning to trait-increasing allele
    if (BETA > 0) {
      EA_I = EA; NEA_I = NEA; EAF_I = EAF
    } else {
      EA_I = NEA; NEA_I = EA; EAF_I = 1 - EAF
    }
    BETA_I = (BETA < 0) ? -BETA : BETA

    # Palindromic SNPs
    PALIND1 = (EA_I == "A" && NEA_I == "T") || (EA_I == "T" && NEA_I == "A")
    PALIND2 = (EA_I == "C" && NEA_I == "G") || (EA_I == "G" && NEA_I == "C")
    PALIND = PALIND1 || PALIND2

    # Alleles in the complementary strand
    comp["A"]="T"
    comp["T"]="A"
    comp["C"]="G"
    comp["G"]="C"
    EA_I_C = comp[EA_I]
    NEA_I_C = comp[NEA_I]

    # Missing EAF
    if (!(EAF_I > 0 && EAF_I < 1)) {
      # Only possible to impute if not palindromic
      if (!PALIND) {
        if ( EA_I == REFA2 && NEA_I == REFA1 ) {
          # Direct match
          EAF_I = REFA2FREQ
        } else if ( EA_I == REFA1 && NEA_I == REFA2 ) {
          # Reverse match
          EAF_I = 1 - REFA2FREQ
        } else if ( EA_I_C == REFA2 && NEA_I_C == REFA1 ) {
          # Complementary match
          EAF_I = REFA2FREQ
        } else if ( EA_I_C == REFA1 && NEA_I_C == REFA2 ) {
          # Reverse-complementary match
          EAF_I = 1 - REFA2FREQ
        } else {
          # No match found, drop
          next
        }
      } else {
        # Palindromic with missing EAF, drop
        next
      }
    }

    # Dropping ambiguous palindromic SNPs (MAF > 0.4)
    if ( PALIND && (EAF_I > 0.4 && EAF_I < 0.6) ) {
      next
    }
    
    # Flipped EAF
    EAF_I_F = 1 - EAF_I

    # Comparing allele frequencies
    # Original
    EAFDIFF1 = EAF_I - REFA2FREQ
    EAFDIFF1 = sqrt(EAFDIFF1*EAFDIFF1)
    # Flipped
    EAFDIFF2 = EAF_I_F - REFA2FREQ
    EAFDIFF2 = sqrt(EAFDIFF2*EAFDIFF2)

    # Harmonizing - tolerating differences in EAF up to 10%
    if ( EA_I == REFA2 && NEA_I == REFA1 && EAFDIFF1 < 0.1 ) {
      # Direct match
      HARMON = 1
    } else if ( EA_I == REFA1 && NEA_I == REFA2 && EAFDIFF2 < 0.1 ) {
      # Reverse match
      HARMON = -1
    } else if ( EA_I_C == REFA2 && NEA_I_C == REFA1 && EAFDIFF1 < 0.1 ) {
      # Complementary match
      HARMON = 1
    } else if ( EA_I_C == REFA1 && NEA_I_C == REFA2 && EAFDIFF2 < 0.1 ) {
      # Reverse complementary match
      HARMON = -1
    } else { 
      # Not possible to harmonize, drop
      next 
    }
    
    # After harmonization, taking alleles from reference panel
    if ( HARMON == 1 ) {
      EA_H = REFA2; NEA_H = REFA1
    } else {
      EA_H = REFA1; NEA_H = REFA2
    }

    # Counting occurrences of chromosome:position
    CHRPOS = $1":"$2
    count[CHRPOS]++

    # Output
    print count[CHRPOS], CHROM, POS, EA_H, NEA_H, EAF_I, BETA_I, SE, PVAL, RSID
    
  }' $1 |\
  awk '$1 == 1' |\
  cut -f2- 
}

# BasalISR
TRAIT=BasalISR
echo "Harmonizing $TRAIT..."
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
{if (($1"\t"$2) in a) print $0,a[$1"\t"$2]}' $REFPANEL - |\
HarmonFx >> $RES

# Glucose
TRAIT=Glucose
echo "Harmonizing $TRAIT..."
SUMSTAT=GCST90002232_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{if (($1"\t"$2) in a) print $1,$2,$3,$4,$5,$6,$7,$8,a[$1"\t"$2]}' $REFPANEL - |\
HarmonFx >> $RES

# GlucoseSens
TRAIT=GlucoseSens
echo "Harmonizing $TRAIT..."
SUMSTAT=ChTEgN
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
unzip -c $SUMDIR/BETACELL/${SUMSTAT} meta1xinsdG30.filthetmafn.rsid.selectedcolumns |\
tail -n+4 |\
sed 's/:/\t/g' |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{if (($1"\t"$2) in a) print $1,$2,toupper($3),toupper($4),$5,$9,$10,$11,a[$1"\t"$2]}' $REFPANEL - |\
HarmonFx >> $RES

# HbA1c
TRAIT=HbA1c
echo "Harmonizing $TRAIT..."
SUMSTAT=GCST90002244_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{if (($1"\t"$2) in a) print $1,$2,$3,$4,$5,$6,$7,$8,a[$1"\t"$2]}' $REFPANEL - |\
HarmonFx >> $RES

# HDL
TRAIT=HDL
echo "Harmonizing $TRAIT..."
SUMSTAT=with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_HDL_INV_ALL_with_N_1.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{if (($2"\t"$3) in a) print $2,$3,$5,$4,$8,$14,$15,$16,a[$2"\t"$3]}' $REFPANEL - |\
HarmonFx >> $RES

# HOMA-IR
TRAIT=HOMA_IR
echo "Harmonizing $TRAIT..."
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
{if (($1"\t"$2) in a) print $0,a[$1"\t"$2]}' $REFPANEL - |\
HarmonFx >> $RES

# Insulin
TRAIT=Insulin
echo "Harmonizing $TRAIT..."
SUMSTAT=GCST90002238_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{if (($1"\t"$2) in a) print $1,$2,$3,$4,$5,$6,$7,$8,a[$1"\t"$2]}' $REFPANEL - |\
HarmonFx >> $RES

# Liver fat
TRAIT=LiverFat
echo "Harmonizing $TRAIT..."
SUMSTAT=GCST90016673_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"}\
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{if (($3"\t"$4) in a) print $3,$4,$5,$6,$7,$8,$9,$2,a[$3"\t"$4]}' $REFPANEL - |\
HarmonFx >> $RES

# Pancreatic fat
TRAIT=PancFat
echo "Harmonizing $TRAIT..."
SUMSTAT=GCST90016675_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"}\
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{if (($3"\t"$4) in a) print $3,$4,$5,$6,$7,$8,$9,$2,a[$3"\t"$4]}' $REFPANEL - |\
HarmonFx >> $RES

# SAT
TRAIT=SAT
echo "Harmonizing $TRAIT..."
SUMSTAT=GCST90016672_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{if (($3"\t"$4) in a) print $3,$4,$5,$6,$7,$8,$9,$2,a[$3"\t"$4]}' $REFPANEL - |\
HarmonFx >> $RES

# TG
TRAIT=TG
echo "Harmonizing $TRAIT..."
SUMSTAT=with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_logTG_INV_ALL_with_N_1.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{if (($2"\t"$3) in a) print $2,$3,$5,$4,$8,$14,$15,$16,a[$2"\t"$3]}' $REFPANEL - |\
HarmonFx >> $RES

# TwoGlucose
TRAIT=TwoGlucose
echo "Harmonizing $TRAIT..."
SUMSTAT=GCST90002227_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{if (($1"\t"$2) in a) print $1,$2,$3,$4,$5,$6,$7,$8,a[$1"\t"$2]}' $REFPANEL - |\
HarmonFx >> $RES

# TwoInsulin
TRAIT=TwoInsulin
echo "Harmonizing $TRAIT..."
SUMSTAT=MAGIC_postchallengeIR_IFC_noBMI_ALL.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
awk 'BEGIN{OFS="\t"}\
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{if (($1"\t"$2) in a) print $1,$2,$3,$4,$7,$5,$6,$8,a[$1"\t"$2]}' $REFPANEL ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
HarmonFx >> $RES

# VAT
TRAIT=VAT
echo "Harmonizing $TRAIT..."
SUMSTAT=GCST90016671_buildGRCh37.tsv.gz
RES=${RESDIR}/${TRAIT}_hrmnzd.tsv
echo -e $HEADER > $RES
gunzip -c ${SUMDIR}/${TRAIT}/${SUMSTAT} |\
awk 'BEGIN{OFS="\t"} \
(FNR==NR){a[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6;next} \
{if (($3"\t"$4) in a) print $3,$4,$5,$6,$7,$8,$9,$2,a[$3"\t"$4]}' $REFPANEL - |\
HarmonFx >> $RES

# Proteins
PQTLPANEL=${DIR}/data/protquery_ukbppp.tsv

tail -n+2 $PQTLPANEL |\
cut -f1,20 |\
sed 's/.tar//g' |\
while read -r GENE GFOLDER
do
   echo "Harmonizing $GENE..."
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
   {if (($1"\t"$3) in a) print $1,$3,$5,$4,$6,$10,$11,10**(-$13),a[$1"\t"$3]}' \
   $REFPANEL - |\
   HarmonFx >> $RES
done