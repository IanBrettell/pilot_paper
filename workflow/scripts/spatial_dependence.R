# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)

# Get variables

## Debug
#IN = "/hps/nobackup/birney/users/ian/pilot/hmm_out/0.08/dist_angle/15.csv"
#N_STATES = 15

## True
IN = snakemake@input[[1]]
OUT_OF = snakemake@output[["open_field"]]
OUT_NO = snakemake@output[["novel_object"]]
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


############################
# Plot OF
############################

FONT_SIZE = 10

ASSAY = "open field"

polar_of = df %>% 
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


spatial_of = df %>% 
  dplyr::filter(assay == ASSAY) %>% 
  dplyr::slice_sample(n = 1e5) %>% 
  ggplot() +
  geom_point(aes(x, y, colour = state_recode),
             alpha = 0.2, size = 0.2) +
  facet_wrap(facets = vars(state_recode),
             nrow = N_ROWS) +
  scale_colour_viridis_c() +
  cowplot::theme_cowplot(font_size = FONT_SIZE) +
  theme(aspect.ratio = 1) +
  ggtitle("spatial distribution") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = "none") +
  xlab("x coordinate") +
  ylab("y coordinate")

of_final = cowplot::ggdraw() +
  cowplot::draw_plot(polar_of,x = 0, y = 0, width = 0.5, height = 1) +
  cowplot::draw_plot(spatial_of,x = 0.5, y = 0, width = 0.5, height = 1) 

ggsave(OUT_OF,
       of_final,
       device = "png",
       width = 16,
       height = 16,
       units = "in",
       dpi = 400)

############################
# Plot NO
############################

FONT_SIZE = 10

ASSAY = "novel object"

polar_no = df %>% 
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


spatial_no = df %>% 
  dplyr::filter(assay == ASSAY) %>% 
  dplyr::slice_sample(n = 1e5) %>% 
  ggplot() +
  geom_point(aes(x, y, colour = state_recode),
             alpha = 0.2, size = 0.2) +
  facet_wrap(facets = vars(state_recode),
             nrow = N_ROWS) +
  scale_colour_viridis_c() +
  cowplot::theme_cowplot(font_size = FONT_SIZE) +
  theme(aspect.ratio = 1) +
  ggtitle("spatial distribution") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = "none") +
  xlab("x coordinate") +
  ylab("y coordinate")

no_final = cowplot::ggdraw() +
  cowplot::draw_plot(polar_no,x = 0, y = 0, width = 0.5, height = 1) +
  cowplot::draw_plot(spatial_no,x = 0.5, y = 0, width = 0.5, height = 1) 

ggsave(OUT_NO,
       no_final,
       device = "png",
       width = 16,
       height = 16,
       units = "in",
       dpi = 400)

