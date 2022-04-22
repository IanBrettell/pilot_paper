# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)

# Get variables

## Debug
IN_FILES = list("/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/0.5/dist_angle/15.rds",
                "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/0.5/dist_angle/20.rds")

## True
IN_FILES = snakemake@input
OUT_PNG = snakemake@output[["png"]]
OUT_PDF = snakemake@output[["pdf"]]

# Read in files

in_list = purrr::map(IN_FILES, readr::read_csv)
