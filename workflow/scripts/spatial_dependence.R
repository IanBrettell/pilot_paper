# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)

# Get variables

## Debug
IN = "/hps/nobackup/birney/users/ian/pilot/hmm_out/0.08/dist_angle/15.csv"
N_STATES = 15

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

# add `line` %>% 
df = df %>% 
  dplyr::mutate(line = dplyr::case_when(fish == "ref" ~ ref_fish,
                                        fish == "test" ~ test_fish)) %>% 
  # recode line
  dplyr::mutate(line = dplyr::recode(line, !!!line_vec)) %>% 
  # order line
  dplyr::mutate(line = factor(line, levels = line_vec)) %>% 
  # recode and order test fish
  dplyr::mutate(test_fish = dplyr::recode(test_fish, !!!line_vec),
                test_fish = factor(test_fish, levels = line_vec))

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
# Plot OF
############################

ASSAY = "open field"
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
  ggtitle(ASSAY) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = "none") +
  xlab("x coordinate") +
  ylab("y coordinate")

############################
# Plot NO
############################

ASSAY = "novel object"
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
  ggtitle(ASSAY) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = "none") +
  xlab("x coordinate") +
  ylab("y coordinate")


############################
# Plot density DGE
############################

# Create density function

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

dens_dge_df = df %>% 
  #dplyr::slice_sample(n = 1e5) %>% 
  # filter out iCab references when paired with a different line
  dplyr::filter(!(fish == "ref" & test_fish != "iCab")) %>% 
  # group by line
  dplyr::group_by(assay, line) %>% 
  # get densities
  dplyr::mutate(density = get_density(x, y, n = 30))

# Run Kruskal-Wallis test

kruskal.test(dens_df_dge$density, dens_df_dge$line)

dens_dge = dens_dge_df %>% 
  # take only states 1:4
  #dplyr::filter(state_recode %in% c(1:4)) %>% 
  # Plot
  ggplot() +
  geom_point(aes(x, y, colour = density),
             alpha = 0.1, size = 0.2) +
  #coord_polar() +
  facet_grid(cols = vars(line), rows = vars(assay)) +
  colorspace::scale_color_continuous_sequential(palette = "Mako", rev = F) +
  #scale_colour_viridis_c(option = "rocket") +
  cowplot::theme_cowplot() +
  theme(aspect.ratio = 1,
        strip.background = element_blank()) +
  xlab("x") +
  ylab("y")

ggsave(here::here("spatial_dens_dge.png"),
       dens_dge,
       device = "png",
       width = 11,
       height = 6,
       units = "in",
       dpi = 400)

############################
# Plot density SGE
############################

dens_sge_df = df %>% 
  # filter out iCab references when paired with a different line
  dplyr::filter(line == "iCab") %>% 
  # group by line
  dplyr::group_by(assay, test_fish) %>% 
  # get densities
  dplyr::mutate(density = get_density(x, y, n = 30)) 

# Run Kruskal Wallis

kruskal.test(dens_sge_df$density, dens_sge_df$test_fish)

# Plot

dens_sge = dens_sge_df %>% 
  # take only states 1:4
  dplyr::filter(state_recode %in% c(1:4)) %>% 
  # Plot
  ggplot() +
  geom_point(aes(x, y, colour = density),
             alpha = 0.1, size = 0.2) +
  #coord_polar() +
  facet_grid(cols = vars(test_fish), rows = vars(assay)) +
  #colorspace::scale_color_continuous_sequential(palette = "Mako", rev = T) +
  scale_colour_viridis_c(option = "rocket") +
  cowplot::theme_cowplot() +
  theme(aspect.ratio = 1,
        strip.background = element_blank()) +
  xlab("x") +
  ylab("y")

ggsave(here::here("spatial_dens_sge.png"),
       dens_sge,
       device = "png",
       width = 11,
       height = 6,
       units = "in",
       dpi = 400)

############################
# Compose final
############################

final = cowplot::plot_grid(dens_dge,
                           dens_sge,
                           nrow = 2,
                           labels = c("A", "B"))

ggsave(OUT,
       final,
       device = "png",
       width = 15,
       height = 12,
       units = "in",
       dpi = 400)

