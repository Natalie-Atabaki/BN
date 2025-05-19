library(readr)
library(dplyr)
library(tidyr)
library(purrr)

DIR <- "/Users/da1078co/Documents/PhD/Projects/BN_Naeimeh"
HITSDIR <- paste0(DIR, "/data/GWAS_Hits")
filelist <- list.files(HITSDIR)
DFCOLTYPES <- "ddccddddcccd"
COMPBASE <- c(
  "A" = "t", 
  "T" = "a", 
  "C" = "g",
  "G" = "c"
)
PGSDIR <- paste0(DIR, "/data/PGS_NoClump")
tibble(
  trait = gsub("_hits\\.tsv", "", filelist),
  file = paste(HITSDIR, filelist, sep = "/")
) %>%
  mutate(
    DAT = map(file, read_tsv, col_types = DFCOLTYPES),
    DAT = map(DAT, arrange, CHROM, POS),
    DAT = map(
      DAT, 
      mutate,
      REFA1_COMP = toupper(stringr::str_replace_all(REFA1, COMPBASE)),
      REFA2_COMP = toupper(stringr::str_replace_all(REFA2, COMPBASE)),
      REFA2FREQ_COMP = 1 - REFA2FREQ,
      MATCH_ORIG = (EA == REFA2 & NEA == REFA1) | (EA == REFA2_COMP & NEA == REFA1_COMP),
      MATCH_FLIP = (EA == REFA1 & NEA == REFA2) | (EA == REFA1_COMP & NEA == REFA2_COMP),
      REFA2FREQ = ifelse(MATCH_ORIG, REFA2FREQ, REFA2FREQ_COMP),
      EAF = coalesce(EAF, REFA2FREQ),
      EA_I = ifelse(BETA > 0, EA, NEA),
      NEA_I = ifelse(BETA > 0, NEA, EA),
      EAF_I = ifelse(BETA > 0, EAF, 1 - EAF)
    ),
    DAT = map(DAT, filter, MATCH_ORIG | MATCH_FLIP),
    DAT = map(
      DAT,
      transmute,
      CHROM, POS, EA = EA_I, NEA = NEA_I, BETA = abs(BETA), SE, PVAL, RSID
    )
  ) %>%
  select(-file) %>%
  unnest(DAT) %>%
  arrange(CHROM, POS, RSID, PVAL) %>%
  distinct(CHROM, POS, RSID, .keep_all = TRUE) %>%
  write_tsv(paste0(PGSDIR, "/NOCLUMP_PGS.tsv"))