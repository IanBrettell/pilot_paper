# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)
library(wesanderson)

# Variables

## Debug

IN = here::here("config/tracking_success.csv")


# Get palette 

line_vec = c("iCab", "HdrR", "HNI", "Kaga", "HO5")
names(line_vec) = c("icab", "hdr", "hni", "kaga", "ho5")
new_pal = grDevices::colorRampPalette(wesanderson::wes_palette("GrandBudapest1"))
pal = new_pal(5)
names(pal) = line_vec

# Read in file

df = readr::read_csv(IN) %>% 
  tidyr::separate(col = "sample",
                  into = c("date", "time", "ref", "test", "tank_side", "quadrant"),
                  sep = "_") %>% 
  dplyr::mutate(test = dplyr::recode(test, !!!line_vec),
                test = factor(test, levels = line_vec)) 


# Plot

df %>% 
  ggplot() +
  geom_histogram(aes(prop_success, fill = test),
                 binwidth = 0.002) +
  cowplot::theme_cowplot() +
  scale_fill_manual(values = pal) +
  xlab("proportion of frames tracked") +
  labs(fill = "test fi")
  
