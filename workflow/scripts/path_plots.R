# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)
library(wesanderson)

# Get variables

## Debug
IN = "/hps/nobackup/birney/users/ian/pilot/hmm_out/0.08/dist_angle/15.csv"
DIMS = here::here("config/split_video_dims.csv")
SAMPLE = "20190615_1305_icab_ho5_R"

## True


###################
# Set parameters
###################

MAX_SECONDS = 180

# Get `darker()` and `lighter()` functions

devtools::source_gist("c5015ee666cdf8d9f7e25fa3c8063c99")

# Get ref and test

split_samp = stringr::str_split(SAMPLE, pattern = "_", simplify = T)

ref = split_samp[3]
test = split_samp[4]

if (ref == test) {
  line_vec = c("iCab_ref", "iCab_test")
  names(line_vec) = c("ref", "test")
  pal = c("#F1BB7B", darker("#F1BB7B", amount = 70))
  names(pal) = line_vec
} else {
  line_vec = c("iCab", "HdrR", "HNI", "Kaga", "HO5")
  names(line_vec) = c("icab", "hdr", "hni", "kaga", "ho5")
  new_pal = grDevices::colorRampPalette(wesanderson::wes_palette("GrandBudapest1"))
  pal = new_pal(5)
  names(pal) = line_vec
}

###################
# Read in files
###################

dims = readr::read_csv(DIMS,
                       col_names = "") %>% 
  # get target rows
  dplyr::filter(sample == SAMPLE)

df = readr::read_csv(IN) %>% 
  # create `sample` variable
  tidyr::unite(col = "sample",
               date, time, ref_fish, test_fish, tank_side,
               sep = "_",
               remove = F) %>% 
  # filter for target sample
  dplyr::filter(sample == SAMPLE)

# Add colour depending on iCab controls or not

if (ref == test){
  df = df %>% 
    dplyr::mutate(line = dplyr::recode(fish, !!!line_vec)) %>% 
    dplyr::mutate(line = factor(line, levels = line_vec))
} else {
  df = df %>% 
    dplyr::mutate(line = dplyr::case_when(fish == "ref" ~ ref_fish,
                                          fish == "test" ~ test_fish)) %>% 
    # recode
    dplyr::mutate(line = dplyr::recode(line, !!!line_vec)) %>% 
    # factorise to order
    dplyr::mutate(line = factor(line, levels = line_vec))
}

# Arrange by assay, frame

df = df %>% 
  dplyr::arrange(assay, line, frame)

###################
# Plot
###################

df %>% 
  dplyr::filter(assay == "open_field" & quadrant == "q1" & seconds > 0 & seconds <= 180) %>% 
  ggplot() +
  geom_path(aes(x, y, colour = line)) +
  scale_colour_manual(values = pal) +
  coord_fixed() +
  cowplot::theme_cowplot() +
  scale_y_reverse()
  
