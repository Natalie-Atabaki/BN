# Functions for MR

suppressMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(survey)
})

# Harmonise GWAS data
harmonfx <- function(joindf) {
  # Complementary base dictionary
  compbase <- c("A" = "t", "T" = "a", "C" = "g", "G" = "c")
  joindf |>
    # Making sure needed columns exist
    select(
      CHROM, POS, RSID,
      EA_EXP, NEA_EXP, EAF_EXP, BETA_EXP, SE_EXP, PVAL_EXP,
      EA_OUT, NEA_OUT, EAF_OUT, BETA_OUT, SE_OUT, PVAL_OUT
    ) |>
    # Removing ambiguous palindromic SNPs with intermediate EAF
    mutate(
      PALIND1 = EA_EXP %in% c("A", "T") & NEA_EXP %in% c("A", "T"),
      PALIND2 = EA_EXP %in% c("C", "G") & NEA_EXP %in% c("C", "G"),
      PALIND = PALIND1 | PALIND2
    ) |>
    filter(
      !(
        PALIND & (
          (EAF_EXP > 0.4 & EAF_EXP < 0.6) |
            (EAF_OUT > 0.4 & EAF_OUT < 0.6)
        )
      )
    ) |>
    mutate(
      # Create complementary alleles
      EA_COMP = stringr::str_replace_all(EA_EXP, compbase) |> toupper(),
      NEA_COMP = stringr::str_replace_all(NEA_EXP, compbase) |> toupper(),
      # Flipped EAF for outcome
      EAF_OUTF = 1 - EAF_OUT,
      HARMON = case_when(
        # Direct match
        EA_EXP == EA_OUT & NEA_EXP == NEA_OUT &
          abs(EAF_EXP - EAF_OUT) < 0.2 ~ 1,
        # Reverse match
        EA_EXP == NEA_OUT & NEA_EXP == EA_OUT &
          abs(EAF_EXP - EAF_OUTF) < 0.2 ~ -1,
        # Direct match, but palindromic
        EA_EXP == EA_OUT & NEA_EXP == NEA_OUT &
          abs(EAF_EXP - EAF_OUTF) < 0.2 & PALIND ~ -1,
        # Reverse match, but palindromic
        EA_EXP == NEA_OUT & NEA_EXP == EA_OUT &
          abs(EAF_EXP - EAF_OUT) < 0.2 & PALIND ~ 1,
        # Match on complement, not palindromic
        EA_COMP == EA_OUT & NEA_COMP == NEA_OUT &
          abs(EAF_EXP - EAF_OUT) < 0.2 & !PALIND ~ 1,
        # Match on reverse complement, not palindromic
        EA_COMP == NEA_OUT & NEA_COMP == EA_OUT &
          abs(EAF_EXP - EAF_OUTF) < 0.2 & !PALIND ~ -1
      ),
      # Harmonizing outcome effect and EAF
      BETA_OUT = BETA_OUT * HARMON,
      EAF_OUT = ifelse(HARMON == 1, EAF_OUT, EAF_OUTF)
    ) |>
    # Keeping only harmonized SNPs with similar EAF
    filter(
      !is.na(HARMON),
      abs(EAF_EXP - EAF_OUT) < 0.2
    ) |>
    # Final table
    transmute(
      CHROM, POS, RSID, EA = EA_EXP, NEA = NEA_EXP,
      EAF_EXP,
      # Capping very small betas and SEs to avoid numerical issues
      BETA_EXP = sign(BETA_EXP) * pmax(abs(BETA_EXP), 1e-22),
      SE_EXP = pmax(SE_EXP, 1e-22),
      PVAL_EXP,
      EAF_OUT,
      # Capping very small betas and SEs to avoid numerical issues
      BETA_OUT = sign(BETA_OUT) * pmax(abs(BETA_OUT), 1e-22),
      SE_OUT = pmax(SE_OUT, 1e-22),
      PVAL_OUT
    )
}

# Clumping
clumpfx <- function(joindf, plinkbin, refpanel) {
  # Temporary file with SNPs to clump
  temploc <- tempfile()
  snpdat <- transmute(joindf, SNP = RSID, P = PVAL_EXP)
  write_tsv(snpdat, temploc)
  # Running PLINK clumping
  system(
    paste(
      plinkbin,
      "--bfile", refpanel,
      "--clump", temploc,
      "--clump-p1", 5e-8,
      "--clump-r2", 0.01,
      "--clump-kb", 250,
      "--silent",
      "--out", temploc
    )
  )
  # Reading clumped SNPs
  clumped <- paste(
    "awk '{print $3}'", 
    paste(temploc, "clumped", sep = ".")
  ) |>
    pipe() |>
    readLines()
  # Keeping only non-empty entries
  clumped <- clumped[nzchar(clumped)]
  # Cleaning up temporary files
  unlink(paste0(temploc, "*"))
  # Returning clumped SNPs
  joindf |>
    filter(RSID %in% clumped)
}

# IVW MR
ivwmrfx <- function(df, plinkbin, refpanel) {
  # Clumping
  clumpdf <- clumpfx(df, plinkbin, refpanel)
  # IVW and Egger regression if more than 1 SNP
  if (nrow(clumpdf) > 1) {
    # Weighted regression - IVW
    wts <- 1 / clumpdf$SE_OUT^2
    ivwreg <- lm(BETA_OUT ~ -1 + BETA_EXP, data = clumpdf, weights = wts)
    ivwreg <- summary(ivwreg)
    # Weighted regression - Egger
    eggerreg <- lm(BETA_OUT ~ BETA_EXP, data = clumpdf, weights = wts)
    eggerreg <- summary(eggerreg)
    # Collecting results
    resdf <- tibble(
      # IVW
      BETA_IVW = ivwreg$coef["BETA_EXP", "Estimate"],
      SE_IVW = ivwreg$coef["BETA_EXP", "Std. Error"],
      # Egger
      BETA_Egger = eggerreg$coef["BETA_EXP", "Estimate"],
      SE_Egger = eggerreg$coef["BETA_EXP", "Std. Error"],
      B0_Egger = eggerreg$coef["(Intercept)", "Estimate"],
      B0SE_Egger = eggerreg$coef["(Intercept)", "Std. Error"],
    ) |>
      mutate(
        # Correcting SEs for overdispersion
        SE_IVW = SE_IVW / min(1, ivwreg$sigma),
        SE_Egger = SE_Egger / min(1, eggerreg$sigma),
        B0SE_Egger = B0SE_Egger / min(1, eggerreg$sigma)
      )
    # Adding p-values and heterogeneity stats
    resdf <- mutate(
      resdf,
      PVAL_IVW = 2 * pnorm(abs(BETA_IVW / SE_IVW), lower.tail = FALSE),
      QDF_IVW = nrow(clumpdf) - 1,
      QSTAT_IVW = (ivwreg$sigma^2) * QDF_IVW,
      QPVAL_IVW = pchisq(QSTAT_IVW, QDF_IVW, lower.tail = FALSE),
      PVAL_Egger = 2 * pnorm(abs(BETA_Egger / SE_Egger), lower.tail = FALSE),
      PVAL_B0Egger = 2 * pnorm(abs(B0_Egger / B0SE_Egger), lower.tail = FALSE),
      QSTAT_Egger = (eggerreg$sigma^2) * (QDF_IVW - 1),
      QPVAL_Egger = pchisq(QSTAT_Egger, QDF_IVW - 1, lower.tail = FALSE),
    )
    datdf <- clumpdf
    failure <- "None"
  } else {
    # If only one SNP, SMR calculation
    resdf <- mutate(
      clumpdf,
      BETA_SMR = BETA_OUT / BETA_EXP,
      SE_SMR = sqrt(
        (BETA_SMR^2) *
          (((SE_EXP / BETA_EXP)^2) + ((SE_OUT / BETA_OUT)^2))
      ),
      PVAL_SMR = pchisq((BETA_SMR / SE_SMR)^2, df = 1, lower.tail = FALSE)
    )
    datdf <- NULL
    failure <- "SingleValidIns"
  }
  # Returning results
  list(RESDF = resdf, DATDF = datdf, FAILURE = failure)
}

# LD matrix calculation
ldmatfx <- function(df, plinkbin, refpanel) {
  # Temporary file with SNPs to get LD for
  ldfile <- tempfile()
  df |> select(RSID, EA) |> write_tsv(ldfile, col_names = FALSE)
  # Running PLINK to get LD matrix
  system(
    paste(
      plinkbin,
      "--bfile", refpanel,
      "--extract", ldfile,
      "--r", "square",
      "--a2-allele", ldfile, 2, 1,
      "--silent",
      "--out", ldfile
    )
  )
  # Reading LD matrix
  ldmat <- read_table(
    paste(ldfile, "ld", sep = "."),
    col_names = FALSE, 
    col_types = cols(.default = "n")
  )
  ldmat <- as.matrix(ldmat)
  colnames(ldmat) <- rownames(ldmat) <- df$RSID
  unlink(paste0(ldfile, "*"))
  ldmat
}

# R2 to top instrument
r2topfx <- function(df, topsnp, plinkbin, refpanel) {
  # Temporary file with SNPs to get R2 for
  r2file <- tempfile()
  df |> select(RSID) |> write_tsv(r2file, col_names = FALSE)
  # Running PLINK to get R2
  system(
    paste(
      plinkbin,
      "--bfile", refpanel,
      "--extract", r2file,
      "--r2",
      "--ld-snp", topsnp, # Reference SNP
      "--ld-window-kb", 1000,
      "--ld-window", 99999999,
      "--ld-window-r2", 0,
      "--silent",
      "--out", r2file
    )
  )
  # Reading R2 table
  r2df <-  read_table(
    paste(r2file, "ld", sep = "."), 
    col_types = '---nncn-',
    skip = 1,
    col_names = c('CHROM', 'POS', 'RSID', 'R2')
  )
  unlink(paste0(r2file, "*"))
  r2df
}

# SMR-HEIDI calculation
smrheidifx <- function(df, minpos, maxpos, plinkbin, refpanel) {
  # Selecting SNPs in the cis region
  cisdf <- filter(df, POS > minpos, POS < maxpos)
  # Top SNP
  topsnpdf <- slice_max(cisdf, BETA_EXP / SE_EXP, n = 1, with_ties = FALSE)
  # SMR calculation
  beta_smr <- with(topsnpdf, BETA_OUT / BETA_EXP)
  se_smr <- with(
    topsnpdf,
    sqrt((beta_smr^2) * (((SE_EXP/BETA_EXP)^2) + ((SE_OUT/BETA_OUT)^2)))
  )
  pval_smr <- pchisq((beta_smr/se_smr)^2, df = 1, lower.tail = FALSE)
  # HEIDI calculation
  # Calculating R2 to top SNP
  r2tab <- r2topfx(df, topsnpdf$RSID, plinkbin, refpanel)
  # Merging R2 info
  df <- inner_join(df, r2tab, by = c("CHROM", "POS", "RSID"))
  # Selecting instruments for HEIDI test
  ins_thres <- pchisq(q = 10, df = 1, lower.tail = FALSE)
  insdf <- df |>
    filter(
      # Top SNP always included
      RSID == topsnpdf$RSID |
        # Other valid instruments in LD with top SNP, only up to 0.9
        (R2 > 0.05 & R2 < 0.9 & PVAL_EXP < ins_thres)
    )
  # HEIDI test only if more than 3 SNPs
  if (nrow(insdf) > 3) {
    # Arranging list of instruments
    insdf <- insdf |>
      arrange(
        # Top SNP first
        desc(RSID == topsnpdf$RSID),
        # Then by R2 to top SNP
        desc(R2),
        # Then by strength of instrument
        desc(BETA_EXP / SE_EXP)
      )
    # Limiting to top SNP + 20 SNPs in LD
    insdf <- slice(insdf, 1:21)
    # Calculating LD matrix
    ldmat <- ldmatfx(insdf, plinkbin, refpanel)
    # Estimated causal effect per instrument
    beta_xy <- with(insdf, BETA_OUT / BETA_EXP)
    # Z-scores of instrument-exposure associations
    zeta_gx <- with(insdf, BETA_EXP / SE_EXP)
    # Covariance of exposure Z-scores
    zeta_gxmat <- zeta_gx %*% t(zeta_gx)
    # Covariance of causal effect estimates
    beta_xy_mat <- beta_xy %*% t(beta_xy)
    # Outer product of standard errors of outcome effects
    se_gymat <- with(insdf, SE_OUT %*% t(SE_OUT))
    # Outer product of exposure effects
    beta_gxmat <- with(insdf, BETA_EXP %*% t(BETA_EXP))
    # Covariance matrix of estimated exposure-outcome effects
    cov_betaxy <- (
      # Uncertainty in outcome effects
      (ldmat * se_gymat / beta_gxmat) +
        # Uncertainty in causal estimates
        (ldmat * beta_xy_mat / zeta_gxmat)
    )
    # Difference between each instrument and top SNP
    diff <- beta_xy - beta_smr
    # Variance of the differences
    var_diff <- cov_betaxy - cov_betaxy[1, ] + (se_smr^2)
    var_diff <- t(t(var_diff) - cov_betaxy[1, ])
    # Removing top SNP from estimated differences
    diff <- diff[-1]
    var_diff <- var_diff[-1, -1]
    # Chi-squared statistic of each difference
    chistat <- (diff^2) / diag(var_diff)
    # Correlation matrix of differences
    corr_diff <- cov2cor(var_diff)
    # Eigenvalues of the correlation matrix
    lambda <- eigen(
      corr_diff, only.values = TRUE,
      symmetric = TRUE # This is a correlation matrix
    )$values
    # HEIDI p-value using the sum of chi-squared statistics
    pval_heidi <- pchisqsum(
      sum(chistat), df = rep(1, length(lambda)),
      a = lambda, method = "saddlepoint", lower.tail = FALSE
    )
    nsnp_heidi <- nrow(insdf) - 1
    fail <- "None"
  } else {
    # Not enough instruments for HEIDI test
    pval_heidi <- NaN
    nsnp_heidi <- NaN
    fail <- "FewValidInsHEIDI"
  }
  # Final results
  resdf <- mutate(
    topsnpdf,
    beta_smr = beta_smr, se_smr = se_smr,
    pval_smr = pval_smr, pval_heidi = pval_heidi,
    nsnp_heidi = nsnp_heidi
  )
  list(resdf = resdf, datdf = df, fail = fail)
}

# Global MR function
mrfx <- function(
  exposure,
  exploc,
  outcome,
  outloc,
  cis = FALSE, chromloc = NA, minpos = NA, maxpos = NA,
  plinkbin = "~/miniforge3/pkgs/plink-1.90b6.21-h2413b67_5/bin/plink",
  refpanel = "~/Documents/Data/1KG/EUR"
) {

  # Check input
  if (cis) {
    stopifnot(!is.na(chromloc), !is.na(minpos), !is.na(maxpos))
  }

  # Reading exposure data
  expdf <- exploc |>
    file.path(paste0(exposure, "_hrmnzd.tsv")) |>
    read_tsv(show_col_types = FALSE) |>
    rename_with(
      function(x) {
        paste(x, "EXP", sep = "_")
      },
      -c(CHROM, POS, RSID)
    ) |>
    # Only valid instruments
    filter(PVAL_EXP < 5e-8)
  # If exposure instruments should be in cis region
  if (cis) {
    expdf <- expdf |>
      filter(
        CHROM == chromloc,
        POS > minpos, POS < maxpos
      )
  }

  # If enough valid exposure instruments
  if (nrow(expdf) > 0) {
    # Proceeding with reading outcome data
    outdf <- outloc |>
      file.path(paste0(outcome, "_hrmnzd.tsv")) |>
      read_tsv(show_col_types = FALSE) |>
      rename_with(
        function(x) {
          paste(x, "OUT", sep = "_")
        },
        -c(CHROM, POS, RSID)
      )
    # Merging exposure and outcome data
    joindf <- inner_join(expdf, outdf, by = c("CHROM", "POS", "RSID"))
    if (nrow(joindf) > 0) {
      # Harmonizing
      joindf <- harmonfx(joindf)
      if (nrow(joindf) > 0) {
        # Running MR
        if (cis) {
          res <- smrheidifx(joindf, minpos, maxpos, plinkbin, refpanel)
        } else {
          res <- ivwmrfx(joindf, plinkbin, refpanel)
        }
      } else {
        # No valid instruments after harmonization
        res <- list(RESDF = NULL, DATDF = NULL, FAILURE = "NoValidInsHarmon")
      }
    } else {
      # No valid instruments after merging exposure and outcome data
      res <- list(RESDF = NULL, DATDF = NULL, FAILURE = "NoMatchForValidIns")
    }
  } else {
    # Not enough valid exposure instruments, returning empty result
    res <- list(RESDF = NULL, DATDF = NULL, FAILURE = "NoValidIns")
  }
  # Returning results
  res
}
