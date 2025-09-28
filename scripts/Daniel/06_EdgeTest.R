setwd('~/')
source('06_UTILS.R')
library(purrr)
edge_table_totest <- read_tsv(here('data/edge_table_totest.tsv'), show_col_types = FALSE)
edge_table_totest <- edge_table_totest |>
  mutate(
    MR_DONE = pmap_lgl(
      list(NODE1, CIS, NODE2, CHROM, start, end),
      function(EXPOSURE, isCIS, OUTCOME, CHROMLOC, STARTPOS, ENDPOS){
        message('Analysing effect of ', EXPOSURE, ' on ', OUTCOME, '...', appendLF = FALSE)
        RESFILE <- paste0(here('data/EdgeTests'), '/', EXPOSURE, '_', OUTCOME, '.rds')
        if(!file.exists(RESFILE)){
          RES <- mrfx(
            EXPOSURE = EXPOSURE, OUTCOME = OUTCOME, 
            CIS = isCIS, CHROMLOC = CHROMLOC, MINPOS = STARTPOS - 1e6, MAXPOS = ENDPOS + 1e6,
            PLINKBIN = '/Users/da1078co/miniforge3/pkgs/plink-1.90b6.21-h2413b67_5/bin/plink',
            REFPANEL = '/Users/da1078co/Documents/Data/1KG/EUR'
          )
          saveRDS(RES, file = RESFILE)
          message('Done.')
        } else {
          message('Link already tested, skipping.')
        }
        return(TRUE)
      }
    )
  )
