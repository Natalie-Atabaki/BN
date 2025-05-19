# Functions for MR

suppressMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(survey)
})

# Harmonise GWAS data
harmonfx <- function(JOINDF){
  COMPBASE <- c("A" = "t", "T" = "a", "C" = "g", "G" = "c")
  JOINDF <- mutate(
    JOINDF,
    EA_COMP = toupper(stringr::str_replace_all(EA_EXP, COMPBASE)),
    NEA_COMP = toupper(stringr::str_replace_all(NEA_EXP, COMPBASE)),
    PALIND1 = EA_EXP %in% c('A', 'T') & NEA_EXP %in% c('A', 'T'),
    PALIND2 = EA_EXP %in% c('C', 'G') & NEA_EXP %in% c('C', 'G'),
    PALIND = PALIND1 | PALIND2,
    EAF_OUTF = 1 - EAF_OUT,
    HARMON = case_when(
      EA_EXP == EA_OUT & NEA_EXP == NEA_OUT & abs(EAF_EXP - EAF_OUT) < 0.2 ~ 1,
      EA_EXP == EA_OUT & NEA_EXP == NEA_OUT & PALIND & abs(EAF_EXP - EAF_OUTF) < 0.2 ~ -1,
      EA_EXP == NEA_OUT & NEA_EXP == EA_OUT & abs(EAF_EXP - EAF_OUTF) < 0.2 ~ -1,
      EA_EXP == NEA_OUT & NEA_EXP == EA_OUT & PALIND & abs(EAF_EXP - EAF_OUT) < 0.2 ~ 1,
      EA_COMP == EA_OUT & NEA_COMP == NEA_OUT & !PALIND & abs(EAF_EXP - EAF_OUT) < 0.2 ~ 1,
      EA_COMP == NEA_OUT & NEA_COMP == EA_OUT & !PALIND & abs(EAF_EXP - EAF_OUTF) < 0.2 ~ -1
    ),
    BETA_OUT = BETA_OUT * HARMON,
    EAF_OUT = ((HARMON == 1) * EAF_OUT) + ((HARMON == -1) * (EAF_OUTF))
  )
  JOINDF <- filter(JOINDF, !is.na(HARMON), abs(EAF_EXP - EAF_OUT) < 0.2)
  JOINDF <- transmute(
    JOINDF, 
    CHROM, POS, RSID, EA = EA_EXP, NEA = NEA_EXP, 
    EAF_EXP, across(c(BETA_EXP, SE_EXP), ~sign(.x) * pmax(abs(.x), 1e-22)), PVAL_EXP,
    EAF_OUT, across(c(BETA_OUT, SE_OUT), ~sign(.x) * pmax(abs(.x), 1e-22)), PVAL_OUT
  )
  return(JOINDF)
}

# Clumping
clumpfx <- function(JOINDF, PLINKBIN, REFPANEL){
  TEMPLOC <- here('data/TEMPFILE')
	SNPDAT <- transmute(JOINDF, SNP = RSID, P = PVAL_EXP)
	write_tsv(SNPDAT, TEMPLOC)
	system(
    paste(
      PLINKBIN,
	    "--bfile", REFPANEL,
	    "--clump", TEMPLOC,
      "--clump-p1", 5e-8,
      "--clump-r2", 0.01,
      "--clump-kb", 250,
      "--silent",
      "--out", TEMPLOC
    )
  )
	CLUMPED <- readLines(
    pipe(paste("awk '{print $3}'", paste(TEMPLOC, "clumped", sep = ".")))
  )
  CLUMPED <- CLUMPED[nzchar(CLUMPED)]
  unlink(paste0(TEMPLOC, "*"))
  RESDF <- filter(JOINDF, RSID %in% CLUMPED)
  return(RESDF)
}

# IVW MR
IVWMRfx <- function(DF, PLINKBIN, REFPANEL){
  CLUMPDF <- clumpfx(DF, PLINKBIN, REFPANEL)
  if(nrow(CLUMPDF) > 1){
    WTS <- 1 / CLUMPDF$SE_OUT^2
    IVWREG <- summary(lm(BETA_OUT ~ -1 + BETA_EXP, data = CLUMPDF, weights = WTS))
    EggerREG <- summary(lm(BETA_OUT ~ BETA_EXP, data = CLUMPDF, weights = WTS))
    RESDF <- tibble(
      BETA_IVW = IVWREG$coef["BETA_EXP","Estimate"],
      SE_IVW = IVWREG$coef["BETA_EXP","Std. Error"] / min(1, IVWREG$sigma),
      BETA_Egger = EggerREG$coef["BETA_EXP","Estimate"],
      SE_Egger = EggerREG$coef["BETA_EXP","Std. Error"] / min(1, EggerREG$sigma),
      B0_Egger = EggerREG$coef["(Intercept)","Estimate"],
      B0SE_Egger = EggerREG$coef["(Intercept)","Std. Error"] / min(1, EggerREG$sigma),
    )
    RESDF <- mutate(
      RESDF,
      PVAL_IVW = 2 * pnorm(abs(BETA_IVW / SE_IVW), lower.tail=FALSE),
      QDF_IVW = nrow(CLUMPDF) - 1,
      QSTAT_IVW = (IVWREG$sigma^2) * QDF_IVW,
      QPVAL_IVW = pchisq(QSTAT_IVW, QDF_IVW, lower.tail=FALSE),
      PVAL_Egger = 2 * pnorm(abs(BETA_Egger / SE_Egger), lower.tail=FALSE),
      PVAL_B0Egger = 2 * pnorm(abs(B0_Egger / B0SE_Egger), lower.tail=FALSE),
      QSTAT_Egger = (EggerREG$sigma^2) * (QDF_IVW - 1),
      QPVAL_Egger = pchisq(QSTAT_Egger, QDF_IVW - 1, lower.tail=FALSE),
    )
    return(list(RESDF = RESDF, DATDF = CLUMPDF, FAILURE = 'None'))
  } else {
    RESDF <- mutate(
      CLUMPDF,
      BETA_SMR = BETA_OUT / BETA_EXP,
      SE_SMR = sqrt((BETA_SMR^2) * (((SE_EXP/BETA_EXP)^2) + ((SE_OUT/BETA_OUT)^2))),
      PVAL_SMR = pchisq((BETA_SMR/SE_SMR)^2, df = 1, lower.tail = FALSE)
    )
    return(list(RESDF = RESDF, DATDF = NULL, FAILURE = 'SingleValidIns'))
  }
}

# LD matrix calculation
LDMATFx <- function(DF, PLINKBIN, REFPANEL){
  LDFILE <- here('data/LDFILE')
  write_tsv(select(DF, RSID, EA), LDFILE)
  system(
    paste(
      PLINKBIN,
      "--bfile", REFPANEL,
      "--extract", LDFILE,
      "--r", "square",
      "--a2-allele", LDFILE, 2, 1,
      "--silent",
      "--out", LDFILE
    )
  )
  LDMAT <- read_table(
    paste(LDFILE, "ld", sep = "."), 
    col_names = FALSE, 
    col_types = cols(.default = "n")
  )
  LDMAT <- as.matrix(LDMAT)
  colnames(LDMAT) <- rownames(LDMAT) <- DF$RSID
  unlink(paste0(LDFILE, "*"))
  return(LDMAT)
}

# R2 to top instrument
R2TopFx <- function(DF, TOPSNP, PLINKBIN, REFPANEL){
  R2FILE <- here('data/R2FILE')
  write_tsv(select(DF, RSID), R2FILE)
  system(
    paste(
      PLINKBIN,
      "--bfile", REFPANEL,
      "--extract", R2FILE,
      "--r2",
      "--ld-snp", TOPSNP,
      "--ld-window-kb", 1000,
      "--ld-window", 99999999, 
      "--ld-window-r2", 0,
      "--silent",
      "--out", R2FILE
    )
  )
  R2DF <-  read_table(
    paste(R2FILE, "ld", sep = "."), 
    col_types = '---nncn-',
    skip = 1,
    col_names = c('CHROM', 'POS', 'RSID', 'R2')
  )
  unlink(paste0(R2FILE, "*"))
  return(R2DF)
}

# SMR-HEIDI calculation
SMRHEIDIfx <- function(DF, MINPOS, MAXPOS, PLINKBIN, REFPANEL){
  CISDF <- filter(DF, POS > MINPOS, POS < MAXPOS)
  TOPSNPDF <- slice_max(CISDF, BETA_EXP/SE_EXP, n = 1, with_ties = FALSE)
  BETA_SMR <- with(TOPSNPDF, BETA_OUT / BETA_EXP)
  SE_SMR <- with(
    TOPSNPDF,
    sqrt((BETA_SMR^2) * (((SE_EXP/BETA_EXP)^2) + ((SE_OUT/BETA_OUT)^2)))
  )
  PVAL_SMR <- pchisq((BETA_SMR/SE_SMR)^2, df = 1, lower.tail = FALSE)
  R2TAB <- R2TopFx(DF, TOPSNPDF$RSID, PLINKBIN, REFPANEL)
  DF <- inner_join(DF, R2TAB, by = c('CHROM', 'POS', 'RSID'))
  INS_THRES <- pchisq(q = 10, df = 1, lower.tail = FALSE)
  INSDF <- filter(DF, RSID == TOPSNPDF$RSID | (R2 > 0.05 & R2 < 0.9 & PVAL_EXP < INS_THRES))
  if(nrow(INSDF) > 3){
    INSDF <- arrange(INSDF, desc(RSID == TOPSNPDF$RSID), desc(R2), desc(BETA_EXP / SE_EXP))
    INSDF <- slice(INSDF, 1:21)
    LDMAT <- LDMATFx(INSDF, PLINKBIN, REFPANEL)
    BETA_XY <- with(INSDF, BETA_OUT / BETA_EXP)
    ZETA_GX <- with(INSDF, BETA_EXP / SE_EXP)
    ZETA_GXMat <- ZETA_GX %*% t(ZETA_GX)
    BETA_XYMat <- BETA_XY %*% t(BETA_XY)
    SE_GYMat <- with(INSDF, SE_OUT %*% t(SE_OUT))
    BETA_GXMat <- with(INSDF, BETA_EXP %*% t(BETA_EXP))
    COV_BETAXY <- ( LDMAT * SE_GYMat / BETA_GXMat ) + ( LDMAT * BETA_XYMat / ZETA_GXMat )
    DIFF <- BETA_XY - BETA_SMR
    VAR_DIFF <- COV_BETAXY - COV_BETAXY[1,] + (SE_SMR^2)
    VAR_DIFF <- t(t(VAR_DIFF) - COV_BETAXY[1,])
    DIFF <- DIFF[-1]
    VAR_DIFF <- VAR_DIFF[-1,-1]
    CHISTAT <- (DIFF^2) / diag(VAR_DIFF)
    m <- nrow(INSDF) - 1
    CORR_DIFF <- diag(m)
    for( i in 1 : (m-1) ) {
      for( j in (i+1) : m ) {
        CORR_DIFF[i,j] = CORR_DIFF[j,i] = (
          VAR_DIFF[i,j] / sqrt(VAR_DIFF[i,i] * VAR_DIFF[j,j])
        )
      }  
    }
    LAMBDA <- eigen(CORR_DIFF, symmetric = TRUE, only.values = TRUE)$values
    PVAL_HEIDI <- pchisqsum(
      sum(CHISTAT), df = rep(1, length(LAMBDA)), 
      a = LAMBDA, method = "saddlepoint", lower.tail = FALSE
    )
    NSNP_HEIDI <- m
    FAIL <- 'None'
  } else {
    PVAL_HEIDI <- NaN
    NSNP_HEIDI <- NaN
    FAIL <- 'FewValidInsHEIDI'
  }
  RESDF <- mutate(
    TOPSNPDF,
    BETA_SMR = BETA_SMR, SE_SMR = SE_SMR, 
    PVAL_SMR = PVAL_SMR, PVAL_HEIDI = PVAL_HEIDI,
    NSNP_HEIDI = NSNP_HEIDI
  )
  return(list(RESDF = RESDF, DATDF = DF, FAILURE = FAIL))
}

# Global MR function
mrfx <- function(
  EXPOSURE, OUTCOME, 
  CIS = FALSE, CHROMLOC = NA, MINPOS = NA, MAXPOS = NA,
  PLINKBIN = '/Users/da1078co/miniforge3/pkgs/plink-1.90b6.21-h2413b67_5/bin/plink',
  REFPANEL = '/Users/da1078co/Documents/Data/1KG/EUR'
){
  
  # Check input
  if(CIS){ stopifnot(!is.na(CHROMLOC), !is.na(MINPOS), !is.na(MAXPOS)) }

  # Exposure data
  EXPFILE <- here(paste0('data/GWAS_Harmonized/', EXPOSURE, '_hrmnzd.tsv'))
  EXPDF <- read_tsv(EXPFILE, show_col_types = FALSE)
  EXPDF <- rename_with(EXPDF, function(x){paste(x, 'EXP', sep = '_')}, -c(CHROM, POS, RSID))
  if(CIS){ 
    EXPDF <- filter(EXPDF, CHROM == CHROMLOC)
    ValidIns <- with(EXPDF, any(POS > MINPOS & POS < MAXPOS & PVAL_EXP < 5e-8))
    if(!ValidIns){ return(list(RESDF = NULL, DATDF = NULL, FAILURE = 'NoValidIns')) }
  } else {
    EXPDF <- filter(EXPDF, PVAL_EXP < 5e-8)
  }
  if(nrow(EXPDF) == 0){ return(list(RESDF = NULL, DATDF = NULL, FAILURE = 'NoValidIns')) }
  
  # Outcome data
  OUTFILE <- here(paste0('data/GWAS_Harmonized/', OUTCOME, '_hrmnzd.tsv'))
  OUTDF <- read_tsv(OUTFILE, show_col_types = FALSE)
  OUTDF <- rename_with(OUTDF, function(x){paste(x, 'OUT', sep = '_')}, -c(CHROM, POS, RSID))
  if(CIS){ OUTDF <- filter(OUTDF, CHROM == CHROMLOC) }
  
  # Joining
  JOINDF <- inner_join(EXPDF, OUTDF, by = c('CHROM', 'POS', 'RSID'))
  if(nrow(JOINDF) == 0){ return(list(RESDF = NULL, DATDF = NULL, FAILURE = 'NoMatchForValidIns')) }

  # Harmonizing
  JOINDF <- harmonfx(JOINDF)
  if(nrow(JOINDF) == 0){ return(list(RESDF = NULL, DATDF = NULL, FAILURE = 'NoValidInsHarmon')) }

  # Running MR
  if(CIS){
    RES <- SMRHEIDIfx(JOINDF, MINPOS, MAXPOS, PLINKBIN, REFPANEL)
  } else {
    RES <- IVWMRfx(JOINDF, PLINKBIN, REFPANEL)
  }

  return(RES)
}