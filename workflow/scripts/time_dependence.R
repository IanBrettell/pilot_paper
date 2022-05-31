# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)
library(ggridges)
library(viridisLite)

# Get variables

## Debug
IN = "/hps/nobackup/birney/users/ian/pilot/hmm_out/0.08/dist_angle/15.csv"
N_STATES = 15

## True
IN = snakemake@input[[1]]
OUT_DGE = snakemake@output[["dge"]]
OUT_SGE = snakemake@output[["sge"]]
N_STATES = snakemake@params[["n_states"]] %>% 
  as.numeric()


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

# Recode `test_fish` and order

df = df %>% 
  dplyr::mutate(test_fish = dplyr::recode(test_fish, !!!line_vec),
                test_fish = factor(test_fish, levels = line_vec))


############################
# Plot polar
############################

#FONT_SIZE = 10
#
#polar_dge = df %>% 
#  # remove iCab when paired with a different test fish
#  dplyr::filter(!(fish == "ref" & test_fish != "icab")) %>% 
#  # select random sample of 1e5 rows
#  dplyr::slice_sample(n = 1e5) %>% 
#  ggplot() +
#  geom_point(aes(angle_recode, log10(distance), colour = state_recode),
#             alpha = 0.3, size = 0.2) +
#  coord_polar() +
#  facet_wrap(~state_recode, nrow = N_ROWS) +
#  scale_x_continuous(labels = c(0, 90, 180, 270),
#                     breaks = c(0, 90, 180, 270)) +
#  scale_color_viridis_c() +
#  guides(colour = "none") +
#  xlab("angle of travel") +
#  ylab(expression(log[10]("distance travelled in pixels"))) +
#  ggtitle("HMM states") +
#  cowplot::theme_cowplot(font_size = FONT_SIZE) +
#  theme(plot.title = element_text(hjust = 0.5)) +
#  facet_wrap(~state_recode, nrow = N_ROWS)
#
#polar_sge = df %>% 
#  # remove iCab when paired with a different test fish
#  dplyr::filter(fish == "ref") %>% 
#  # select random sample of 1e5 rows
#  dplyr::slice_sample(n = 1e5) %>% 
#  ggplot() +
#  geom_point(aes(angle_recode, log10(distance), colour = state_recode),
#             alpha = 0.3, size = 0.2) +
#  coord_polar() +
#  facet_wrap(~state_recode, nrow = N_ROWS) +
#  scale_x_continuous(labels = c(0, 90, 180, 270),
#                     breaks = c(0, 90, 180, 270)) +
#  scale_color_viridis_c(option = "inferno") +
#  guides(colour = "none") +
#  xlab("angle of travel") +
#  ylab(expression(log[10]("distance travelled in pixels"))) +
#  #ggtitle("HMM states") +
#  cowplot::theme_cowplot(font_size = FONT_SIZE) +
#  #theme(plot.title = element_text(hjust = 0.5)) +
#  facet_wrap(~state_recode, nrow = N_ROWS)
#

############################
# Plot ridges
############################

#ridge_plot = df %>% 
#  # remove iCab when paired with a different test fish
#  dplyr::filter(!(fish == "ref" & test_fish != "icab")) %>% 
#  # filter for target assay
#  #dplyr::filter(assay == ASSAY) %>% 
#  dplyr::mutate(state_recode = factor(state_recode, levels = 1:N_STATES)) %>% 
#  ggplot() +
#  ggridges::geom_density_ridges(aes(x = seconds, y = state_recode, fill = state_recode)) +
#  scale_y_discrete(limits = rev) +
#  facet_grid(rows = vars(assay),
#             cols = vars(line)) +
#  scale_fill_viridis_d() +
#  cowplot::theme_cowplot(font_size = FONT_SIZE) +
#  guides(fill = "none") +
#  ylab("HMM state") +
#  scale_x_continuous(breaks = c(0,200,400,600))

# Set significant states

#######################
# Medarkov matrices: DGE
#######################

SEC_INT = 2

dge_tile_df = df %>% 
  # remove iCab ref fishes (because DGE compares test fishes)
  dplyr::filter(!(line == "iCab" & fish == "ref")) %>% 
  #dplyr::slice_sample(n = 1e6) %>% 
  # add `indiv` column
  tidyr::unite(date, time, quadrant, fish,
               col = "indiv",
               sep = "_") %>%
  # get top state per 2 seconds
  dplyr::mutate(seconds_bin = floor(seconds / SEC_INT)) %>% 
  dplyr::group_by(assay, indiv, line, seconds_bin) %>% 
  dplyr::count(state_recode) %>% 
  dplyr::slice_max(order_by = n, n = 1) %>% 
  dplyr::ungroup() %>% 
  # reverse order by `indiv` so that the earliest videos are at the top
  dplyr::arrange(indiv) %>% 
  # convert `seconds_bin` back to seconds
  dplyr::mutate(seconds = seconds_bin * SEC_INT)

# Open field

dge_tile_of = dge_tile_df %>% 
  dplyr::filter(assay == "open field") %>% 
  ggplot() +
  geom_tile(aes(seconds, indiv, fill = state_recode)) + 
  facet_grid(rows = vars(line), cols = vars(assay), scales = "free") +
  scale_fill_viridis_c() +
  scale_y_discrete(limits = rev) +
  guides(fill = "none") +
  ylab("individual fish") +
  cowplot::theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) 

# Novel object

dge_tile_no = dge_tile_df %>% 
  dplyr::filter(assay == "novel object") %>% 
  ggplot() +
  geom_tile(aes(seconds, indiv, fill = state_recode)) + 
  facet_grid(rows = vars(line), cols = vars(assay), scales = "free") +
  scale_fill_viridis_c() +
  scale_y_discrete(limits = rev) +
  guides(fill = "none") +
  ylab("individual fish") +
  cowplot::theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) 

#######################
# Medarkov matrices: SGE
#######################

sge_tile_df = df %>% 
  # take only iCab test fishes 
  dplyr::filter(fish == "ref") %>% 
  # add `indiv` column
  tidyr::unite(date, time, quadrant, fish,
               col = "indiv",
               sep = "_") %>%
  # rename and reorder assay
  dplyr::mutate(assay = stringr::str_replace(assay, "_", " "),
                assay = factor(assay, levels = c("open field", "novel object"))) %>% 
  # get top state per 2 seconds
  dplyr::mutate(seconds_bin = floor(seconds / SEC_INT)) %>% 
  dplyr::group_by(assay, indiv, test_fish, seconds_bin) %>% 
  dplyr::count(state_recode) %>% 
  dplyr::slice_max(order_by = n, n = 1) %>% 
  dplyr::ungroup() %>% 
  # reverse order by `indiv` so that the earliest videos are at the top
  dplyr::arrange(indiv) %>% 
  # convert `seconds_bin` back to seconds
  dplyr::mutate(seconds = seconds_bin * SEC_INT)

# Open field

sge_tile_of = sge_tile_df %>% 
  dplyr::filter(assay == "open field") %>% 
  ggplot() +
  geom_tile(aes(seconds, indiv, fill = state_recode)) + 
  facet_grid(rows = vars(test_fish), cols = vars(assay), scales = "free") +
  scale_fill_viridis_c(option = "inferno") +
  scale_y_discrete(limits = rev) +
  guides(fill = "none") +
  ylab("individual fish") +
  cowplot::theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) 

# Novel object

sge_tile_no = sge_tile_df %>% 
  dplyr::filter(assay == "novel object") %>% 
  ggplot() +
  geom_tile(aes(seconds, indiv, fill = state_recode)) + 
  facet_grid(rows = vars(test_fish), cols = vars(assay), scales = "free") +
  scale_fill_viridis_c(option = "inferno") +
  scale_y_discrete(limits = rev) +
  guides(fill = "none") +
  ylab("individual fish") +
  cowplot::theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) 

##########################
# Time density - DGE
##########################

SIGS_DGE_OF = c(1,2,3,5)
SIGS_DGE_NO = c(1,2,3)

# Take viridis colours for significant states and add grey
pal_dge_of = viridisLite::viridis(n = N_STATES)
pal_dge_of = c(pal_dge_of[SIGS_DGE_OF], "#9da2ab")
names(pal_dge_of) = c(as.character(SIGS_DGE_OF), "other")

pal_dge_no = viridisLite::viridis(n = N_STATES)
pal_dge_no = c(pal_dge_no[SIGS_DGE_NO], "#9da2ab")
names(pal_dge_no) = c(as.character(SIGS_DGE_NO), "other")

ASSAY = "open field"
time_dens_dge_of = df %>% 
  dplyr::filter(assay == ASSAY) %>% 
  # remove iCab when paired with a different test fish
  dplyr::filter(fish == "test") %>% 
  # filter for target assay
  #dplyr::filter(assay == ASSAY) %>% 
  # recode state 
  dplyr::mutate(state_plot_recode = dplyr::case_when(state_recode %in% SIGS_DGE_OF ~ as.character(state_recode),
                                                     TRUE ~ "other"),
                state_plot_recode = factor(state_plot_recode, levels = c(as.character(SIGS_DGE_OF), "other"))) %>% 
  ggplot() +
  geom_density(aes(seconds, after_stat(count), fill = state_plot_recode),
               position = "fill") +
  facet_grid(rows = vars(line),
             cols = vars(assay)) + 
  scale_fill_manual(values = pal_dge_of) +
  cowplot::theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  guides(fill = "none") +
  scale_x_continuous(breaks = c(0,200,400,600)) +
  scale_y_continuous(breaks = c(0,0.5,1))

ASSAY = "novel object"
time_dens_dge_no = df %>% 
  dplyr::filter(assay == ASSAY) %>% 
  # remove iCab when paired with a different test fish
  dplyr::filter(fish == "test") %>% 
  # filter for target assay
  #dplyr::filter(assay == ASSAY) %>% 
  # recode state 
  dplyr::mutate(state_plot_recode = dplyr::case_when(state_recode %in% SIGS_DGE_NO ~ as.character(state_recode),
                                                     TRUE ~ "other"),
                state_plot_recode = factor(state_plot_recode, levels = c(as.character(SIGS_DGE_NO), "other"))) %>% 
  ggplot() +
  geom_density(aes(seconds, after_stat(count), fill = state_plot_recode),
               position = "fill") +
  facet_grid(rows = vars(line),
             cols = vars(assay)) + 
  scale_fill_manual(values = pal_dge_no) +
  cowplot::theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  guides(fill = "none") +
  scale_x_continuous(breaks = c(0,200,400,600)) +
  scale_y_continuous(breaks = c(0,0.5,1))


##########################
# Time density: SGE
##########################

SIGS_SGE_OF = 1:3
SIGS_SGE_NO = 1:3

# Take viridis colours for significant states and add grey
pal_sge_of = viridisLite::viridis(n = N_STATES, option = "inferno")
pal_sge_of = c(pal_sge_of[SIGS_SGE_OF], "#9da2ab")
names(pal_sge_of) = c(as.character(SIGS_SGE_OF), "other")

pal_sge_no = viridisLite::viridis(n = N_STATES, option = "inferno")
pal_sge_no = c(pal_sge_no[SIGS_SGE_NO], "#9da2ab")
names(pal_sge_no) = c(as.character(SIGS_SGE_NO), "other")

ASSAY = "open field"
time_dens_sge_of = df %>% 
  dplyr::filter(assay == ASSAY) %>% 
  # remove iCab when paired with a different test fish
  dplyr::filter(fish == "ref") %>% 
  # filter for target assay
  #dplyr::filter(assay == ASSAY) %>% 
  # recode state 
  dplyr::mutate(state_plot_recode = dplyr::case_when(state_recode %in% SIGS_SGE_OF ~ as.character(state_recode),
                                                     TRUE ~ "other"),
                state_plot_recode = factor(state_plot_recode, levels = c(as.character(SIGS_SGE_OF), "other"))) %>% 
  ggplot() +
  geom_density(aes(seconds, after_stat(count), fill = state_plot_recode),
               position = "fill") +
  facet_grid(rows = vars(test_fish),
             cols = vars(assay)) + 
  scale_fill_manual(values = pal_sge_of) +
  cowplot::theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  guides(fill = "none") +
  scale_x_continuous(breaks = c(0,200,400,600)) +
  scale_y_continuous(breaks = c(0,0.5,1))

ASSAY = "novel object"
time_dens_sge_no = df %>% 
  dplyr::filter(assay == ASSAY) %>% 
  # remove iCab when paired with a different test fish
  dplyr::filter(fish == "ref") %>% 
  # filter for target assay
  #dplyr::filter(assay == ASSAY) %>% 
  # recode state 
  dplyr::mutate(state_plot_recode = dplyr::case_when(state_recode %in% SIGS_SGE_NO ~ as.character(state_recode),
                                                     TRUE ~ "other"),
                state_plot_recode = factor(state_plot_recode, levels = c(as.character(SIGS_SGE_NO), "other"))) %>% 
  ggplot() +
  geom_density(aes(seconds, after_stat(count), fill = state_plot_recode),
               position = "fill") +
  facet_grid(rows = vars(test_fish),
             cols = vars(assay)) + 
  scale_fill_manual(values = pal_sge_no) +
  cowplot::theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  guides(fill = "none") +
  scale_x_continuous(breaks = c(0,200,400,600)) +
  scale_y_continuous(breaks = c(0,0.5,1))

##########################
# Spatial density: DGE
##########################

# Create density function

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


# OF

sdens_dge_df = df %>% 
  #dplyr::slice_sample(n = 1e5) %>% 
  # filter out iCab references when paired with a different line
  dplyr::filter(fish == "test") %>% 
  # group by line
  dplyr::group_by(assay, line) %>% 
  # get densities
  dplyr::mutate(density = get_density(x, y, n = 30)) %>% 
  dplyr::ungroup()

  # Plot
sdens_dge_of = sdens_dge_df %>% 
  # take only states 1:4
  dplyr::filter(state_recode %in% SIGS_DGE_OF & assay == "open field") %>% 
  ggplot() +
  geom_point(aes(x, y, colour = density),
             alpha = 0.1, size = 0.2) +
  facet_grid(cols = vars(assay), rows = vars(line)) +
  colorspace::scale_color_continuous_sequential(palette = "Mako", rev = F) +
  #scale_colour_viridis_c(option = "rocket") +
  cowplot::theme_cowplot(font_size = 12) +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  scale_x_continuous(breaks = c(0,150,300)) +
  scale_y_continuous(breaks = c(0,150,300)) +
  xlab("x") +
  ylab("y")

# NO

sdens_dge_no = sdens_dge_df %>% 
  dplyr::filter(state_recode %in% SIGS_DGE_NO & assay == "novel object") %>% 
  # Plot
  ggplot() +
  geom_point(aes(x, y, colour = density),
             alpha = 0.1, size = 0.2) +
  #coord_polar() +
  facet_grid(cols = vars(assay), rows = vars(line)) +
  colorspace::scale_color_continuous_sequential(palette = "Mako", rev = F) +
  #scale_colour_viridis_c(option = "rocket") +
  cowplot::theme_cowplot(font_size = 12) +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  scale_x_continuous(breaks = c(0,150,300)) +
  scale_y_continuous(breaks = c(0,150,300)) +
  xlab("x") +
  ylab("y")

##########################
# Spatial density: SGE
##########################

# OF

sdens_sge_df = df %>% 
  #dplyr::slice_sample(n = 1e5) %>% 
  # filter out iCab references when paired with a different line
  dplyr::filter(fish == "ref") %>% 
  # group by line
  dplyr::group_by(assay, test_fish) %>% 
  # get densities
  dplyr::mutate(density = get_density(x, y, n = 30)) %>% 
  dplyr::ungroup()

# Plot
sdens_sge_of = sdens_sge_df %>% 
  # take only states 1:4
  dplyr::filter(state_recode %in% SIGS_SGE_OF & assay == "open field") %>% 
  ggplot() +
  geom_point(aes(x, y, colour = density),
             alpha = 0.1, size = 0.2) +
  facet_grid(cols = vars(assay), rows = vars(test_fish)) +
  #colorspace::scale_color_continuous_sequential(palette = "Mako", rev = F) +
  scale_colour_viridis_c(option = "rocket") +
  cowplot::theme_cowplot(font_size = 12) +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  scale_x_continuous(breaks = c(0,150,300)) +
  scale_y_continuous(breaks = c(0,150,300)) +
  xlab("x") +
  ylab("y")

# NO

sdens_sge_no = sdens_sge_df %>% 
  dplyr::filter(state_recode %in% SIGS_SGE_NO & assay == "novel object") %>% 
  # Plot
  ggplot() +
  geom_point(aes(x, y, colour = density),
             alpha = 0.1, size = 0.2) +
  #coord_polar() +
  facet_grid(cols = vars(assay), rows = vars(test_fish)) +
  #colorspace::scale_color_continuous_sequential(palette = "Mako", rev = F) +
  scale_colour_viridis_c(option = "rocket") +
  cowplot::theme_cowplot(font_size = 12) +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  scale_x_continuous(breaks = c(0,150,300)) +
  scale_y_continuous(breaks = c(0,150,300)) +
  xlab("x") +
  ylab("y")


############################
# Combine and save
############################

#final = cowplot::ggdraw() +
#  cowplot::draw_plot(polar_plot,
#                     x = 0, y = 0,
#                     width = 0.4, height = 1) +
#  cowplot::draw_plot(ridge_plot,
#                     x = 0.4, y = 0,
#                     width = 0.6, height = 1) 

#top = cowplot::plot_grid(polar_dge, time_dens_dge,
#                         rel_widths = c(0.333, 0.666),
#                         axis = c("bt"),
#                         align = c("hv"))
#
#bottom = cowplot::plot_grid(polar_sge, time_dens_sge,
#                            rel_widths = c(0.333, 0.666),
#                            axis = c("bt"),
#                            align = c("hv"))
#
#final = cowplot::plot_grid(top, bottom,
#                           rel_widths = c(1,1),
#                           rel_heights = c(0.5, 0.5),
#                           nrow = 2)


final_dge = cowplot::plot_grid(dge_tile_of +
                                 theme(strip.background.y = element_blank(),
                                       strip.text.y = element_blank(),
                                       axis.title.y = element_text(vjust=-5)),
                               time_dens_dge_of,
                               sdens_dge_of,
                               dge_tile_no +
                                 theme(strip.background.y = element_blank(),
                                       strip.text.y = element_blank(),
                                       axis.title.y = element_text(vjust=-5)),
                               time_dens_dge_no,
                               sdens_dge_no,
                               nrow = 2, ncol = 3,
                               rel_widths = c(1,1,0.6,1,1,0.6),
                               align = "hv",
                               labels = c('A', 'B', 'C', 'D', 'E', 'F'))


final_sge = cowplot::plot_grid(sge_tile_of +
                                 theme(strip.background.y = element_blank(),
                                       strip.text.y = element_blank(),
                                       axis.title.y = element_text(vjust=-5)),
                               time_dens_sge_of,
                               sdens_sge_of,
                               sge_tile_no +
                                 theme(strip.background.y = element_blank(),
                                       strip.text.y = element_blank(),
                                       axis.title.y = element_text(vjust=-5)),
                               time_dens_sge_no,
                               sdens_sge_no,
                               nrow = 2, ncol = 3,
                               rel_widths = c(1,1,0.6,1,1,0.6),
                               align = "hv",
                               labels = c('A', 'B', 'C', 'D', 'E', 'F'))


ggsave(OUT_DGE,
       final_dge,
       device = "png",
       width = 11.5,
       height = 12,
       units = "in",
       dpi = 400)

ggsave(OUT_SGE,
       final_sge,
       device = "png",
       width = 11.5,
       height = 12,
       units = "in",
       dpi = 400)

