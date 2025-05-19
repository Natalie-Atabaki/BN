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
PLINKBIN <- "/Users/da1078co/miniforge3/pkgs/plink-1.90b6.21-h2413b67_5/bin/plink"
REFPANEL <- "/Users/da1078co/RefData/1KG/EUR"
PGSDIR <- paste0(DIR, "/data/PGS_Scores")
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
      REFA2FREQ_FLIP = 1 - REFA2FREQ,
      MATCH_ORIG = (EA == REFA2 & NEA == REFA1) | (EA == REFA2_COMP & NEA == REFA1_COMP),
      MATCH_FLIP = (EA == REFA1 & NEA == REFA2) | (EA == REFA1_COMP & NEA == REFA2_COMP),
      REFA2FREQ = ifelse(MATCH_ORIG, REFA2FREQ, REFA2FREQ_FLIP),
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
    ),
    DAT = map(  
      DAT,
      function(DF){
        TEMPLOC <- paste(PGSDIR, "TEMPFILE", sep = "/")
	SNPDAT <- transmute(DF, SNP = RSID, P = PVAL)
	write_tsv(SNPDAT, TEMPLOC)
	system(
	  paste(
	    PLINKBIN,
	    "--bfile", REFPANEL,
	    "--clump", TEMPLOC,
            "--clump-p1", 1e-6,
            "--clump-r2", 0.01,
            "--clump-kb", 250,
            "--out", TEMPLOC
          )
        )
	CLUMPED <- readLines(
	  pipe(paste("awk '{print $3}'", paste(TEMPLOC, "clumped", sep = ".")))
        )
        CLUMPED <- CLUMPED[nzchar(CLUMPED)]
        unlink(paste0(TEMPLOC, "*"))
        return(filter(DF, RSID %in% CLUMPED))
      }
    ),
    DAT = map(DAT, arrange, CHROM, POS),
    DAT = map2_lgl(
      DAT, trait,
      function(DF, TRAIT){
        write_tsv(DF, paste0(PGSDIR, "/", TRAIT, "_PGS.tsv"))
        return(TRUE)
      }
    )
  )