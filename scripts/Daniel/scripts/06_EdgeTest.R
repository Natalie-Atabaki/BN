# Working directory
setwd("/Users/da1078co/Documents/Lund/PhD/Projects/BN")
# Libraries and functions
source("scripts/Daniel/scripts/00_Functions.R")
library(purrr)
# Load edges to test
edge_table_totest <- read_tsv(
  file.path("data", "edge_table_totest.tsv"),
  show_col_types = FALSE
)
# Running analyses
edge_table_totest <- edge_table_totest |>
  mutate(
    MR_DONE = pmap_lgl(
      list(node1, node2, atype, CHROM, start, end),
      function(
        exposure, outcome,
        analysistype, chromloc, startpos, endpos
      ) {
        message(
          "Analysing effect of ", exposure, " on ",
          outcome, "...", appendLF = FALSE
        )
        resfile <- file.path(
          "data", "EdgeTests",
          paste0(exposure, "_", outcome, ".rds")
        )
        res <- mrfx(
          exposure, exploc = file.path("data", "GWAS_Harmonized"),
          outcome, outloc = file.path("data", "GWAS_Harmonized"),
          analysistype,
          chromloc,
          minpos = startpos - 1e6, maxpos = endpos + 1e6,
          plinkbin = file.path(
            "/Users/da1078co/miniforge3",
            "pkgs/plink-1.90b6.21-h2413b67_5/bin/plink"
          ),
          refpanel = "/Users/da1078co/Documents/Data/1KG/EUR"
        )
        write_rds(res, file = resfile)
        message("Done.")
        TRUE
      }
    )
  )
