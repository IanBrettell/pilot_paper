# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(wesanderson)
#library(ggbeeswarm)
library(rstatix)
#library(ggpubr)
library(cowplot)

# Get variables

## Debug
IN = "/hps/nobackup/birney/users/ian/pilot/hmm_out/0.08/dist_angle/15.csv"
N_STATES = 15
VARIABLES = "distance and angle of travel"
INTERVAL = 0.08

# Create line recode vector
line_vec = c("iCab", "HdrR", "HNI", "Kaga", "HO5")
names(line_vec) = c("icab", "hdr", "hni", "kaga", "ho5")


SEC_INT = 2

# Read in file

df = readr::read_csv(IN) %>% 
  # recode angle to sit between 0 and 360
  dplyr::mutate(angle_recode = ifelse(angle < 0,
                                      180 + (180 + angle),
                                      angle))

# Recode states by mean distance

rank_df = df %>% 
  dplyr::group_by(state) %>% 
  dplyr::summarise(mean_dist = mean(distance)) %>% 
  # rank
  dplyr::arrange(mean_dist) %>% 
  dplyr::mutate(rank = 1:nrow(.))

recode_vec = rank_df %>% 
  dplyr::pull(rank)
names(recode_vec) = rank_df %>% 
  dplyr::pull(state)

# Recode `state`

df = df %>% 
  dplyr::mutate(state_recode = dplyr::recode(state, !!!recode_vec))

# Get SGE df

sge_df = df %>% 
  # keep only iCabs
  ## add `line`
  dplyr::mutate(line = dplyr::case_when(fish == "ref" ~ ref_fish,
                                        fish == "test" ~ test_fish)) %>% 
  # recode `test_fish`
  dplyr::mutate(test_fish = dplyr::recode(test_fish, !!!line_vec)) %>% 
  # order `test_fish`
  dplyr::mutate(test_fish = factor(test_fish, levels = line_vec)) %>% 
  # rename and reorder assay
  dplyr::mutate(assay = stringr::str_replace(assay, "_", " "),
                assay = factor(assay, levels = c("open field", "novel object"))) %>% 
  # add `run` column
  tidyr::unite(date, time, quadrant,
               col = "run",
               sep = "_",
               remove = F) %>% 
  # add `indiv` column
  tidyr::unite(date, time, quadrant, fish,
               col = "indiv",
               sep = "_")

#######################
# Co-occupancy boxplot
#######################

# Per 2 seconds
cooc_bin = sge_df %>% 
  # get top state per 2 seconds
  dplyr::mutate(seconds_bin = floor(seconds / SEC_INT)) %>% 
  dplyr::group_by(run, assay, indiv, test_fish, seconds_bin) %>% 
  dplyr::count(state_recode) %>% 
  dplyr::slice_max(order_by = n, n = 1) %>% 
  dplyr::ungroup() 

cooc = sge_df %>% 
  # pivot wider to get cols for ref and test
  tidyr::pivot_wider(id_cols = c("run", "assay", "test_fish", "seconds"),
                     names_from = fish,
                     values_from = state_recode) %>% 
  # group by run and assay
  dplyr::group_by(run, assay, test_fish) %>% 
  summarise(TOTAL_ROWS = n(),
            TOTAL_CONC = sum(ref == test, na.rm = T),
            FREQ_CONC = TOTAL_CONC / TOTAL_ROWS) %>% 
  dplyr::ungroup()

# Plot

cooc %>% 
  ggplot() +
  geom_boxplot(aes(test_fish, FREQ_CONC)) +
  facet_grid(cols = vars(assay))

#######################
# Co-occupancy heatmap
#######################


