# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)
library(ggridges)

# Get variables

## Debug
#IN = "/hps/nobackup/birney/users/ian/pilot/hmm_out/0.08/dist_angle/15.csv"
#N_STATES = 15

## True
IN = snakemake@input[[1]]
OUT = snakemake@output[["fig"]]
N_STATES = snakemake@params[["n_states"]]


#######################
# Read in data
#######################

# Get number of rows (for plotting) based on number of states

if (N_STATES == 15){
  N_ROWS = 5
} else if (N_STATES == 12 | 16) {
  N_ROWS = 4
} else if (N_STATES == 17 | 18){
  N_ROWS = 6
}


# Create line recode vector
line_vec = c("iCab", "HdrR", "HNI", "Kaga", "HO5")
names(line_vec) = c("icab", "hdr", "hni", "kaga", "ho5")

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

# Recode `assay`

df = df %>% 
  dplyr::mutate(assay = stringr::str_replace(assay, "_", " "),
                assay = factor(assay, levels = c("open field", "novel object")))

# Add `line`

df = df %>% 
  dplyr::mutate(line = dplyr::case_when(fish == "ref" ~ ref_fish,
                                        fish == "test" ~ test_fish)) %>% 
  # recode
  dplyr::mutate(line = dplyr::recode(line, !!!line_vec)) %>% 
  # factorise to order
  dplyr::mutate(line = factor(line, levels = line_vec))


############################
# Plot polar
############################

FONT_SIZE = 10

polar_plot = df %>% 
  # select random sample of 1e5 rows
  dplyr::slice_sample(n = 1e5) %>% 
  ggplot() +
  geom_point(aes(angle_recode, log10(distance), colour = state_recode),
             alpha = 0.3, size = 0.2) +
  coord_polar() +
  facet_wrap(~state_recode, nrow = N_ROWS) +
  scale_x_continuous(labels = c(0, 90, 180, 270),
                     breaks = c(0, 90, 180, 270)) +
  scale_color_viridis_c() +
  guides(colour = "none") +
  xlab("angle of travel") +
  ylab(expression(log[10]("distance travelled in pixels"))) +
  ggtitle("distance and angle") +
  cowplot::theme_cowplot(font_size = FONT_SIZE) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~state_recode, nrow = N_ROWS)

############################
# Plot ridges
############################

ridge_plot = df %>% 
  # remove iCab when paired with a different test fish
  dplyr::filter(!(fish == "ref" & ref_fish == "icab" & test_fish != "icab")) %>% 
  # filter for target assay
  #dplyr::filter(assay == ASSAY) %>% 
  dplyr::mutate(state_recode = factor(state_recode, levels = 1:N_STATES)) %>% 
  ggplot() +
  ggridges::geom_density_ridges(aes(x = seconds, y = state_recode, fill = state_recode)) +
  scale_y_discrete(limits = rev) +
  facet_grid(rows = vars(assay),
             cols = vars(line)) +
  scale_fill_viridis_d() +
  cowplot::theme_cowplot(font_size = FONT_SIZE) +
  guides(fill = "none") +
  ylab("HMM state") +
  scale_x_continuous(breaks = c(0,200,400,600))
  
  
############################
# Combine and save
############################

final = cowplot::ggdraw() +
  cowplot::draw_plot(polar_plot,
                     x = 0, y = 0,
                     width = 0.4, height = 1) +
  cowplot::draw_plot(ridge_plot,
                     x = 0.4, y = 0,
                     width = 0.6, height = 1) 

ggsave(OUT,
       final,
       device = "png",
       width = 24,
       height = 16,
       units = "in",
       dpi = 400)


