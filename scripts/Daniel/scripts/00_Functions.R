# Functions for MR

# Libraries
library(readr)
library(dplyr)
library(tidyr)
library(survey)

#' Clumping function using PLINK
#' @param joindf Data frame with SNP exposure and outcome association data.
#' @returns Data frame with data harmonized to the exposure alleles.
#' @details
#' Input data have previously harmonized to a reference panel.
#' @importFrom dplyr mutate case_when filter rename transmute across all_of
harmonfx <- function(joindf) {
  joindf %>%
    mutate(
      EAF_OUTF = 1 - .data$EAF_OUT,
      HARMON = case_when(
        # Direct match
        EA_EXP == EA_OUT & NEA_EXP == NEA_OUT &
          abs(EAF_EXP - EAF_OUT) < 0.1 ~ 1,
        # Reverse match
        EA_EXP == NEA_OUT & NEA_EXP == EA_OUT &
          abs(EAF_EXP - EAF_OUTF) < 0.1 ~ -1
        # Only two options due to prior harmonization
      ),
      # Harmonizing outcome effect and EAF
      BETA_OUT = .data$BETA_OUT * .data$HARMON,
      EAF_OUT = ifelse(.data$HARMON == 1, .data$EAF_OUT, .data$EAF_OUTF)
    ) %>%
    # Keeping only harmonized SNPs with similar EAF
    filter(!is.na(.data$HARMON)) %>%
    rename(EA = .data$EA_EXP, NEA = .data$NEA_EXP) %>%
    # Final table
    transmute(
      across(
        all_of(
          c(
            "CHROM", "POS", "RSID",
            "EA", "NEA", "EAF_EXP", "EAF_OUT",
            "BETA_EXP", "SE_EXP", "PVAL_EXP",
            "BETA_OUT", "SE_OUT", "PVAL_OUT"
          )
        )
      )
    )
}


#' Clumping function using PLINK
#' @param joindf Data frame with SNPs to clump.
#' @param plinkbin Path to PLINK binary.
#' @param refpanel Path to reference panel (PLINK fileset).
#' @param clump_pval Clumping p-value threshold.
#' @param clump_r2 Clumping r2 threshold.
#' @param clump_kb Clumping kb window.
#' @returns Data frame with clumped SNPs.
#' @importFrom dplyr filter transmute
#' @importFrom readr write_tsv
clumpfx <- function(
  joindf,
  plinkbin,
  refpanel,
  clump_pval,
  clump_r2,
  clump_kb
) {
  # Temporary file with SNPs to clump
  temploc <- tempfile()
  snpdat <- transmute(
    joindf,
    SNP = .data$RSID,
    P = .data$PVAL_EXP
  )
  write_tsv(snpdat, temploc)
  # Running PLINK clumping
  system(
    paste(
      plinkbin,
      "--bfile", refpanel,
      "--clump", temploc,
      "--clump-p1", clump_pval,
      "--clump-r2", clump_r2,
      "--clump-kb", clump_kb,
      "--silent",
      "--out", temploc
    )
  )
  # Reading clumped SNPs
  clumped <- paste(
    "awk '{print $3}'",
    paste(temploc, "clumped", sep = ".")
  ) %>%
    pipe() %>%
    readLines()
  # Keeping only non-empty entries
  clumped <- clumped[nzchar(clumped)]
  # Cleaning up temporary files
  unlink(paste0(temploc, "*"))
  # Returning clumped SNPs
  joindf %>%
    filter(.data$RSID %in% clumped)
}

#-----------------#
# IVW-MR function #
#-----------------#

#' Inverse-variance weighted and Egger MR
#' @param df Data frame with harmonized SNP exposure and outcome associations
#' @param plinkbin Path to PLINK binary.
#' @param refpanel Path to reference panel (PLINK fileset).
#' @param clump_pval Clumping p-value threshold.
#' @param clump_r2 Clumping r2 threshold.
#' @param clump_kb Clumping kb window.
#' @returns
#' A list with 3 elements:
#' * `resdf` A tibble with IVW and Egger MR results.
#' * `datdf` A tibble with the clumped SNPs used as instruments.
#' * `failure` A character string indicating if there were errors.
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
ivwmrfx <- function(
  df,
  plinkbin,
  refpanel,
  clump_pval,
  clump_r2,
  clump_kb
) {
  # Clumping
  clumpdf <- clumpfx(
    df, plinkbin, refpanel,
    clump_pval, clump_r2, clump_kb
  )
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
      beta_ivw = ivwreg$coef["BETA_EXP", "Estimate"],
      se_ivw = ivwreg$coef["BETA_EXP", "Std. Error"],
      # Egger
      beta_egger = eggerreg$coef["BETA_EXP", "Estimate"],
      se_egger = eggerreg$coef["BETA_EXP", "Std. Error"],
      b0_egger = eggerreg$coef["(Intercept)", "Estimate"],
      b0se_egger = eggerreg$coef["(Intercept)", "Std. Error"],
    ) %>%
      mutate(
        # Correcting SEs for overdispersion
        se_ivw = .data$se_ivw / min(1, ivwreg$sigma),
        se_egger = .data$se_egger / min(1, eggerreg$sigma),
        b0se_egger = .data$b0se_egger / min(1, eggerreg$sigma)
      )
    # Adding p-values and heterogeneity stats
    resdf <- mutate(
      resdf,
      pval_ivw = 2 * pnorm(
        abs(.data$beta_ivw / .data$se_ivw),
        lower.tail = FALSE
      ),
      qdf_ivw = nrow(clumpdf) - 1,
      qstat_ivw = (ivwreg$sigma^2) * .data$qdf_ivw,
      qpval_ivw = pchisq(.data$qstat_ivw, .data$qdf_ivw, lower.tail = FALSE),
      pval_egger = 2 * pnorm(
        abs(.data$beta_egger / .data$se_egger),
        lower.tail = FALSE
      ),
      pval_b0egger = 2 * pnorm(
        abs(.data$b0_egger / .data$b0se_egger),
        lower.tail = FALSE
      ),
      qstat_egger = (eggerreg$sigma^2) * (.data$qdf_ivw - 1),
      qpval_egger = pchisq(
        .data$qstat_egger, .data$qdf_ivw - 1,
        lower.tail = FALSE
      ),
    )
    datdf <- clumpdf
    failure <- "None"
  } else {
    # If only one SNP, SMR calculation
    resdf <- mutate(
      clumpdf,
      beta_smr = .data$BETA_OUT / .data$BETA_EXP,
      se_smr = sqrt(
        (.data$BETA_SMR^2) * (
          ((.data$SE_EXP / .data$BETA_EXP)^2) +
            ((.data$SE_OUT / .data$BETA_OUT)^2)
        )
      ),
      pval_smr = pchisq(
        (.data$beta_smr / .data$se_smr)^2, df = 1, lower.tail = FALSE
      )
    )
    datdf <- NULL
    failure <- "SingleValidIns"
  }
  # Returning results
  list(resdf = resdf, datdf = datdf, failure = failure)
}

#---------------------#
# SMR-HEIDI functions #
#---------------------#

#' LD matrix calculation
#' @param df Data frame with SNPs.
#' @param plinkbin Path to PLINK binary.
#' @param refpanel Path to reference panel (PLINK fileset).
#' @returns LD matrix.
#' @importFrom readr read_table write_tsv
#' @importFrom dplyr select
ldmatfx <- function(df, plinkbin, refpanel) {
  # Temporary file with SNPs to get LD for
  ldfile <- tempfile()
  df %>%
    select(all_of("RSID", "EA")) %>%
    write_tsv(ldfile, col_names = FALSE)
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

#' R2 to top instrument
#' @param df Data frame with SNPs.
#' @param topsnp RSID of the top SNP.
#' @param plinkbin Path to PLINK binary.
#' @param refpanel Path to reference panel (PLINK fileset).
#' @param kbwindow KB window for R2 calculation.
#' @returns R2 data frame.
#' @importFrom readr read_table write_tsv
#' @importFrom dplyr select
r2topfx <- function(df, topsnp, plinkbin, refpanel, kbwindow) {
  # Temporary file with SNPs to get R2 for
  r2file <- tempfile()
  df %>%
    select(.data$RSID) %>%
    write_tsv(r2file, col_names = FALSE)
  # Running PLINK to get R2
  system(
    paste(
      plinkbin,
      "--bfile", refpanel,
      "--extract", r2file,
      "--r2",
      "--ld-snp", topsnp, # Reference SNP
      "--ld-window-kb", kbwindow,
      "--ld-window", 99999999,
      "--ld-window-r2", 0,
      "--silent",
      "--out", r2file
    )
  )
  # Reading R2 table
  r2df <-  read_table(
    paste(r2file, "ld", sep = "."),
    col_types = "-----cn-",
    skip = 1,
    col_names = c("RSID", "R2")
  )
  unlink(paste0(r2file, "*"))
  r2df
}

#' SMR-HEIDI calculation
#' @param df Data frame with SNP exposure and outcome associations.
#' @param minpos Minimum position for cis region.
#' @param maxpos Maximum position for cis region.
#' @param ins_thres P-value to select SNPs for HEIDI. Default: pchisq=10/df=1.
#' @param plinkbin Path to PLINK binary.
#' @param refpanel Path to reference panel (PLINK fileset).
#' @param kbwindow KB window for R2 calculation.
#' @returns
#' A list with 3 elements:
#' * `resdf` A tibble with SMR and HEIDI results.
#' * `datdf` A tibble with the SNPs used as instruments.
#' * `failure` A character string indicating if there were errors.
#' @importFrom readr read_table write_tsv
#' @importFrom dplyr filter mutate slice_min inner_join arrange desc
#' @importFrom survey pchisqsum
smrheidifx <- function(
  df,
  minpos, maxpos,
  ins_thres = pchisq(q = 10, df = 1, lower.tail = FALSE),
  plinkbin, refpanel,
  kbwindow
) {
  # Selecting SNPs in the cis region
  cisdf <- df %>%
    filter(
      .data$POS > minpos,
      .data$POS < maxpos
    ) %>%
    # Cap extremely small outcome betas to avoid numerical issues
    mutate(
      BETA_OUT = ifelse(.data$BETA_OUT < 0, -1, 1) *
        pmax(abs(.data$BETA_OUT), 1e-10)
    )
  # Top SNP
  topsnpdf <- slice_min(cisdf, .data$PVAL_EXP, n = 1, with_ties = FALSE)
  # SMR calculation
  beta_smr <- with(topsnpdf, BETA_OUT / BETA_EXP)
  se_smr <- with(
    topsnpdf,
    sqrt(
      (beta_smr^2) * (
        ((SE_EXP / BETA_EXP)^2) +
          ((SE_OUT / BETA_OUT)^2)
      )
    )
  )
  pval_smr <- pchisq((beta_smr / se_smr)^2, df = 1, lower.tail = FALSE)
  # HEIDI calculation
  # Calculating R2 to top SNP
  r2tab <- r2topfx(
    cisdf, topsnpdf$RSID,
    plinkbin, refpanel,
    kbwindow
  )
  # Merging R2 info
  cisdf <- cisdf %>%
    inner_join(r2tab, by = "RSID")
  # Selecting instruments for HEIDI test
  insdf <- cisdf %>%
    filter(
      # Top SNP always included
      .data$RSID == topsnpdf$RSID |
        # Other valid instruments in LD with top SNP
        (
          .data$R2 > 0.05 & # Excluding weak LD
            .data$R2 < 0.9 & # Excluding perfect LD
            .data$PVAL_EXP < ins_thres # Sufficiently strong instruments
        )
    )
  # HEIDI test only if more than 3 SNPs
  if (nrow(insdf) > 3) {
    # Arranging list of instruments
    insdf <- insdf %>%
      arrange(
        # Top SNP first
        desc(.data$RSID == topsnpdf$RSID),
        # Then by R2 to top SNP
        desc(.data$R2),
        # Then by strength of instrument
        desc(.data$BETA_EXP / .data$SE_EXP)
      ) %>%
      # Limiting to top SNP + 20 SNPs in LD
      slice(1:21)
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
  list(resdf = resdf, datdf = cisdf, failure = fail)
}

#' Global MR function
#' @param exposure String with exposure GWAS name.
#' @param exploc Parent folder for exposure GWAS harmonized data.
#' @param outcome String with outcome GWAS name.
#' @param outloc Parent folder for outcome GWAS harmonized data.
#' @param atype Analysis type, either "SMRHEIDI" or "IVWMR".
#' @param plinkbin Path to PLINK binary.
#' @param refpanel Path to reference panel (PLINK fileset).
#' @param pval_thres (`atype="IVWMR"`) P-value to select exposure instruments.
#' @param clump_r2 (`atype="IVWMR"`) Clumping r2 threshold.
#' @param clump_kb (`atype="IVWMR"`) Clumping kb window.
#' @param chromloc (`atype="SMRHEIDI"`) Chromosome for cis region.
#' @param minpos (`atype="SMRHEIDI"`) Minimum position for cis region.
#' @param maxpos (`atype="SMRHEIDI"`) Maximum position for cis region.
#' @param kbwindow (`atype="SMRHEIDI"`) KB window to measure R2 to top SNP.
#' @param ins_thres (`atype="SMRHEIDI"`) P-value to select SNPs for HEIDI.
#' @returns
#' A list with 3 elements:
#' * `resdf` A tibble with SMR and HEIDI results.
#' * `datdf` A tibble with the SNPs used as instruments.
#' * `failure` A character string indicating if there were errors.
#' @importFrom readr read_tsv
#' @importFrom dplyr filter inner_join
mrfx <- function(
  exposure,
  exploc,
  outcome,
  outloc,
  atype,
  plinkbin,
  refpanel,
  pval_thres = 5e-8,
  clump_r2 = 0.01,
  clump_kb = 250,
  chromloc = NA, minpos = NA, maxpos = NA,
  kbwindow = 1000,
  ins_thres = pchisq(q = 10, df = 1, lower.tail = FALSE)
) {
  # Reading exposure data
  expdf <- exploc %>%
    file.path(paste0(exposure, "_hrmnzd.tsv")) %>%
    read_tsv(col_types = "nnccnnnnc") %>%
    # Only valid instruments
    filter(.data$PVAL < pval_thres)
  # If exposure instruments should be in cis region
  if (atype == "SMRHEIDI") {
    expdf <- expdf %>%
      filter(
        .data$CHROM == chromloc,
        .data$POS > minpos,
        .data$POS < maxpos
      )
  }
  # If enough valid exposure instruments
  if (nrow(expdf) > 0) {
    # Proceeding to read outcome data
    outdf <- outloc %>%
      file.path(paste0(outcome, "_hrmnzd.tsv")) %>%
      read_tsv(col_types = "nnccnnnnc")
    # Merging exposure and outcome data
    joindf <- inner_join(
      expdf, outdf,
      by = c("CHROM", "POS", "RSID"),
      suffix = c("_EXP", "_OUT")
    )
    if (nrow(joindf) > 0) {
      # Harmonizing
      joindf <- harmonfx(joindf)
      if (nrow(joindf) > 0) {
        # Running MR
        res <- switch(
          atype,
          SMRHEIDI = smrheidifx(
            joindf, minpos, maxpos, ins_thres,
            plinkbin, refpanel, kbwindow
          ),
          IVWMR = ivwmrfx(
            joindf, plinkbin, refpanel,
            pval_thres, clump_r2, clump_kb
          ),
        )
      } else {
        # No valid instruments after harmonization
        res <- list(resdf = NULL, datdf = NULL, failure = "NoValidInsHarmon")
      }
    } else {
      # No valid instruments after merging exposure and outcome data
      res <- list(resdf = NULL, datdf = NULL, failure = "NoMatchForValidIns")
    }
  } else {
    # Not enough valid exposure instruments, returning empty result
    res <- list(resdf = NULL, datdf = NULL, failure = "NoValidIns")
  }
  # Returning results
  res
}
