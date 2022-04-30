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

## True
IN = snakemake@input[[1]]
N_STATES = snakemake@params[["n_states"]] %>% 
  as.numeric()
INTERVAL = snakemake@params[["interval"]] %>% 
  as.numeric()
VARIABLES = "distance and angle of travel"
POLAR_BOX_DGE = snakemake@output[["polar_box_dge"]]
POLAR_BOX_SGE = snakemake@output[["polar_box_sge"]]
POLAR_BOX_DGE_SGE = snakemake@output[["polar_box_dge_sge"]]
TILE_DGE = snakemake@output[["tile_dge"]]
TILE_SGE = snakemake@output[["tile_sge"]]
TILE_DGE_SGE = snakemake@output[["tile_dge_sge"]]

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

# Get figure height

HEIGHT = 2 * N_ROWS
 
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

#######################
# Polar plots
#######################

# Set title
TITLE = paste("N states: ",
              N_STATES,
              "\nVariables: ",
              VARIABLES,
              "\nInterval: ",
              INTERVAL, " seconds",
              sep = "")

polar_dge = df %>% 
  # remove iCab when paired with a different test fish
  dplyr::filter(!(fish == "ref" & ref_fish == "icab" & test_fish != "icab")) %>% 
  # select random sample of 1e5 rows
  dplyr::slice_sample(n = 1e5) %>% 
  # factorise `state_recode`
  #dplyr::mutate(state_recode = factor(state_recode, levels = recode_vec)) %>% 
  ggplot() +
  geom_point(aes(angle_recode, log10(distance), colour = state_recode),
             alpha = 0.3, size = 0.2) +
  coord_polar() +
  facet_wrap(~state_recode, nrow = N_ROWS) +
  theme_dark(base_size = 8) +
  scale_x_continuous(labels = c(0, 90, 180, 270),
                     breaks = c(0, 90, 180, 270)) +
  scale_color_viridis_c() +
  guides(colour = "none") +
  xlab("angle of travel") +
  ylab(expression(log[10]("distance travelled in pixels"))) +
  ggtitle(TITLE)

polar_sge = df %>% 
  # keep only iCabs
  ## add `line`
  dplyr::mutate(line = dplyr::case_when(fish == "ref" ~ ref_fish,
                                        fish == "test" ~ test_fish)) %>% 
  dplyr::filter(line == "icab") %>% 
  # select random sample of 1e5 rows
  dplyr::slice_sample(n = 1e5) %>% 
  # factorise `state_recode`
  #dplyr::mutate(state_recode = factor(state_recode, levels = recode_vec)) %>% 
  ggplot() +
  geom_point(aes(angle_recode, log10(distance), colour = state_recode),
             alpha = 0.3, size = 0.2) +
  coord_polar() +
  facet_wrap(~state_recode, nrow = N_ROWS) +
  theme_dark(base_size = 8) +
  scale_x_continuous(labels = c(0, 90, 180, 270),
                     breaks = c(0, 90, 180, 270)) +
  scale_color_viridis_c(option = "magma") +
  guides(colour = "none") +
  xlab("angle of travel") +
  ylab(expression(log[10]("distance travelled in pixels"))) +
  ggtitle(TITLE)

#######################
# Frequency plot
#######################

# Create palette

new_pal = grDevices::colorRampPalette(wesanderson::wes_palette("GrandBudapest1"))
pal = new_pal(5)
names(pal) = line_vec

# DGE

dge_df = df %>% 
  # remove iCab when paired with a different test fish
  dplyr::filter(!(fish == "ref" & ref_fish == "icab" & test_fish != "icab")) %>% 
  # add `line` %>% 
  dplyr::mutate(line = dplyr::case_when(fish == "ref" ~ ref_fish,
                                        fish == "test" ~ test_fish)) %>% 
  # recode line
  dplyr::mutate(line = dplyr::recode(line, !!!line_vec)) %>% 
  # order line
  dplyr::mutate(line = factor(line, levels = line_vec))

#dge_df %>% 
#  # get frequency of time spent in each state for each line
#  dplyr::group_by(line) %>% 
#  dplyr::count(line, state_recode) %>% 
#  dplyr::add_count(line, wt = n, name = "nn") %>% 
#  dplyr::mutate(state_freq = n / nn) %>% 
#  ggplot() +
#  geom_col(aes(state_recode, state_freq, fill = state_recode)) +
#  facet_grid(rows = vars(line)) +
#  scale_fill_viridis_c() +
#  theme_bw() +
#  scale_x_continuous(breaks = unique(dge_df$state_recode)) +
#  guides(fill = "none") +
#  xlab("HMM state") +
#  ylab("Proportion of time spent in HMM state")

# SGE

sge_df = df %>% 
  # keep only iCabs
  ## add `line`
  dplyr::mutate(line = dplyr::case_when(fish == "ref" ~ ref_fish,
                                        fish == "test" ~ test_fish)) %>% 
  dplyr::filter(line == "icab") %>% 
  # recode `test_fish`
  dplyr::mutate(test_fish = dplyr::recode(test_fish, !!!line_vec)) %>% 
  # order `test_fish`
  dplyr::mutate(test_fish = factor(test_fish, levels = line_vec))

#sge_df %>% 
#  # get frequency of time spent in each state for each line
#  dplyr::group_by(test_fish) %>% 
#  dplyr::count(test_fish, state_recode) %>% 
#  dplyr::add_count(test_fish, wt = n, name = "nn") %>% 
#  dplyr::mutate(state_freq = n / nn) %>% 
#  ggplot() +
#  geom_col(aes(state_recode, state_freq, fill = state_recode)) +
#  facet_grid(rows = vars(test_fish)) +
#  scale_fill_viridis_c(option = "magma") +
#  theme_bw() +
#  scale_x_continuous(breaks = unique(dge_df$state_recode)) +
#  guides(fill = "none") +
#  xlab("HMM state") +
#  ylab("Proportion of time spent in HMM state")

#######################
# Box plots -- DGE
#######################  

# Get proportions of time spent in each state

state_freq_dge = dge_df %>% 
  ## unite columns to get reads for each fish
  tidyr::unite(date, time, quadrant, fish,
               col = "indiv",
               sep = "_") %>% 
  ## count rows per fish per state
  dplyr::count(indiv, assay, line, state_recode) %>% 
  # add total row count per fish
  dplyr::add_count(indiv, assay, line, wt = n, name = "nn") %>% 
  # get proportion of time fish spent in each state
  dplyr::mutate(state_freq = n / nn)
  
# Run Kruskal-Wallis test

kw_dge = state_freq_dge %>% 
  dplyr::group_by(assay, state_recode) %>% 
  rstatix::kruskal_test(state_freq ~ line) %>% 
  rstatix::adjust_pvalue(method = "fdr") %>% 
  dplyr::mutate(p.adj = signif(p.adj, digits = 3)) %>% 
  rstatix::add_significance(p.col = "p.adj") %>% 
  # paste p-value and significance together
  dplyr::mutate(p_final = dplyr::case_when(p.adj.signif == "ns" ~ paste("p =", p.adj),
                                           TRUE ~ paste("p =", p.adj, p.adj.signif)))


# Plot boxplots

## Open field

ASSAY = "open_field"
of_box_dge = state_freq_dge %>% 
  dplyr::filter(assay == ASSAY) %>% 
  ggplot() +
    geom_boxplot(aes(line, state_freq, fill = line), notch = T) +
    #ggbeeswarm::geom_beeswarm(aes(line, state_freq), size = 0.5) +
    #ggpubr::stat_pvalue_manual(data = kw_dge_2, label = "p.adj", remove.bracket = T, y.position = 0.8) +
    geom_text(data = kw_dge %>% 
                dplyr::filter(assay == ASSAY),
              aes(x = "HdrR", y = 0.7, label = p_final),
              size = 3) +
    scale_fill_manual(values = pal) +
    theme_bw() +
    facet_wrap(~state_recode, nrow = N_ROWS) +
    guides(fill = "none") +
    ylab("HMM state frequency") +
    ggtitle(stringr::str_replace(ASSAY, "_", " ")) +
    theme(plot.title = element_text(hjust = 0.5))

## Novel object

ASSAY = "novel_object"
no_box_dge = state_freq_dge %>% 
  dplyr::filter(assay == ASSAY) %>% 
  ggplot() +
  geom_boxplot(aes(line, state_freq, fill = line), notch = T) +
  #ggbeeswarm::geom_beeswarm(aes(line, state_freq), size = 0.5) +
  #ggpubr::stat_pvalue_manual(data = kw_dge_2, label = "p.adj", remove.bracket = T, y.position = 0.8) +
  geom_text(data = kw_dge %>% 
              dplyr::filter(assay == ASSAY),
            aes(x = "HdrR", y = 0.7, label = p_final),
            size = 3) +
  scale_fill_manual(values = pal) +
  theme_bw() +
  facet_wrap(~state_recode, nrow = N_ROWS) +
  guides(fill = "none") +
  ylab("HMM state frequency") +
  ggtitle(stringr::str_replace(ASSAY, "_", " ")) +
  theme(plot.title = element_text(hjust = 0.5))

#######################
# Box plots -- SGE
#######################  

# Get proportions of time spent in each state

state_freq_sge = sge_df %>% 
  ## unite columns to get reads for each fish
  tidyr::unite(date, time, quadrant, fish,
               col = "indiv",
               sep = "_") %>% 
  ## count rows per fish per state
  dplyr::count(indiv, assay, test_fish, state_recode) %>% 
  # add total row count per fish
  dplyr::add_count(indiv, assay, test_fish, wt = n, name = "nn") %>% 
  # get proportion of time fish spent in each state
  dplyr::mutate(state_freq = n / nn)

# Run Kruskal-Wallis test

kw_sge = state_freq_sge %>% 
  dplyr::group_by(assay, state_recode) %>% 
  rstatix::kruskal_test(state_freq ~ test_fish) %>% 
  rstatix::adjust_pvalue(method = "fdr") %>% 
  dplyr::mutate(p.adj = signif(p.adj, digits = 3)) %>% 
  rstatix::add_significance(p.col = "p.adj") %>% 
  # paste p-value and significance together
  dplyr::mutate(p_final = dplyr::case_when(p.adj.signif == "ns" ~ paste("p =", p.adj),
                                           TRUE ~ paste("p =", p.adj, p.adj.signif)))


# Plot boxplots

## Open field

ASSAY = "open_field"
of_box_sge = state_freq_sge %>% 
  dplyr::filter(assay == ASSAY) %>% 
  ggplot() +
  geom_boxplot(aes(test_fish, state_freq, colour = test_fish), notch = T) +
  #ggbeeswarm::geom_beeswarm(aes(line, state_freq), size = 0.5) +
  #ggpubr::stat_pvalue_manual(data = kw_dge_2, label = "p.adj", remove.bracket = T, y.position = 0.8) +
  geom_text(data = kw_sge %>% 
              dplyr::filter(assay == ASSAY),
            aes(x = "HdrR", y = 0.7, label = p_final),
            size = 3) +
  #scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  theme_bw() +
  facet_wrap(~state_recode, nrow = N_ROWS) +
  guides(colour = "none") +
  xlab("line of tank partner") +
  ylab("HMM state frequency") +
  ggtitle(stringr::str_replace(ASSAY, "_", " ")) +
  theme(plot.title = element_text(hjust = 0.5))

## Novel object

ASSAY = "novel_object"
no_box_sge = state_freq_sge %>% 
  dplyr::filter(assay == ASSAY) %>% 
  ggplot() +
  geom_boxplot(aes(test_fish, state_freq, colour = test_fish), notch = T) +
  #ggbeeswarm::geom_beeswarm(aes(line, state_freq), size = 0.5) +
  #ggpubr::stat_pvalue_manual(data = kw_dge_2, label = "p.adj", remove.bracket = T, y.position = 0.8) +
  geom_text(data = kw_sge %>% 
              dplyr::filter(assay == ASSAY),
            aes(x = "HdrR", y = 0.7, label = p_final),
            size = 3) +
  scale_colour_manual(values = pal) +
  theme_bw() +
  facet_wrap(~state_recode, nrow = N_ROWS) +
  guides(colour = "none") +
  xlab("line of tank partner") +
  ylab("HMM state frequency") +
  ggtitle(stringr::str_replace(ASSAY, "_", " ")) +
  theme(plot.title = element_text(hjust = 0.5))


#######################
# Compile final polar + boxplots
#######################

# DGE

FONT_SIZE = 10
final_dge = cowplot::ggdraw() +
  cowplot::draw_plot(polar_dge + 
                       ggtitle("HMM states") + 
                       cowplot::theme_cowplot(font_size = FONT_SIZE) +
                       theme(plot.title = element_text(hjust = 0.5)) +
                       facet_wrap(~state_recode, nrow = N_ROWS),
                     x = 0, y = 0,
                     width = 0.4,
                     height = 1) +
  cowplot::draw_plot(of_box_dge +
                      cowplot::theme_cowplot(font_size = FONT_SIZE) +
                       theme(plot.title = element_text(hjust = 0.5)) +
                      facet_wrap(~state_recode, nrow = N_ROWS),
                     x = 0.4, y = 0,
                     width = 0.3,
                     height = 1) +
  cowplot::draw_plot(no_box_dge +
                       cowplot::theme_cowplot(font_size = FONT_SIZE) +
                       theme(plot.title = element_text(hjust = 0.5)) +
                       facet_wrap(~state_recode, nrow = N_ROWS),
                     x = 0.7, y = 0,
                     width = 0.3,
                     height = 1) +
  cowplot::draw_plot_label(label = c("A", "B", "C"),
                           x = c(0, 0.4, 0.7),
                           y = c(1, 1, 1))

#ggsave(here::here("tmp.png"),
#       final_dge,
#       device = "png",
#       width = 16.8,
#       height = 12,
#       units = "in",
#       dpi = 400)

ggsave(POLAR_BOX_DGE,
       final_dge,
       device = "png",
       width = 16.8,
       height = HEIGHT,
       units = "in",
       dpi = 400)

# SGE

final_sge = cowplot::ggdraw() +
  cowplot::draw_plot(polar_sge + 
                       ggtitle("HMM states") + 
                       cowplot::theme_cowplot(font_size = FONT_SIZE) +
                       theme(plot.title = element_text(hjust = 0.5)) +
                       facet_wrap(~state_recode, nrow = N_ROWS),
                     x = 0, y = 0,
                     width = 0.4,
                     height = 1) +
  cowplot::draw_plot(of_box_sge +
                       cowplot::theme_cowplot(font_size = FONT_SIZE) +
                       theme(plot.title = element_text(hjust = 0.5)) +
                       facet_wrap(~state_recode, nrow = N_ROWS),
                     x = 0.4, y = 0,
                     width = 0.3,
                     height = 1) +
  cowplot::draw_plot(no_box_sge +
                       cowplot::theme_cowplot(font_size = FONT_SIZE) +
                       theme(plot.title = element_text(hjust = 0.5)) +
                       facet_wrap(~state_recode, nrow = N_ROWS),
                     x = 0.7, y = 0,
                     width = 0.3,
                     height = 1) +
  cowplot::draw_plot_label(label = c("A", "B", "C"),
                           x = c(0, 0.4, 0.7),
                           y = c(1, 1, 1))

ggsave(POLAR_BOX_SGE,
       final_sge,
       device = "png",
       width = 16.8,
       height = HEIGHT,
       units = "in",
       dpi = 400)

# Together

final_dge_sge = cowplot::ggdraw() +
  cowplot::draw_plot(polar_dge + 
                       ggtitle("HMM states") + 
                       cowplot::theme_cowplot(font_size = FONT_SIZE) +
                       theme(plot.title = element_text(hjust = 0.5)) +
                       facet_wrap(~state_recode, nrow = N_ROWS),
                     x = 0, y = 0.5,
                     width = 0.333,
                     height = 0.5) +
  cowplot::draw_plot(of_box_dge +
                       cowplot::theme_cowplot(font_size = FONT_SIZE) +
                       theme(plot.title = element_text(hjust = 0.5)) +
                       facet_wrap(~state_recode, nrow = N_ROWS),
                     x = 0.333, y = 0.5,
                     width = 0.333,
                     height = 0.5) +
  cowplot::draw_plot(no_box_dge +
                       cowplot::theme_cowplot(font_size = FONT_SIZE) +
                       theme(plot.title = element_text(hjust = 0.5)) +
                       facet_wrap(~state_recode, nrow = N_ROWS),
                     x = 0.666, y = 0.5,
                     width = 0.333,
                     height = 0.5) +
  cowplot::draw_plot_label(label = c("A", "B", "C"),
                           x = c(0, 0.333, 0.666),
                           y = c(1, 1, 1)) +
  cowplot::draw_plot(polar_sge + 
                       ggtitle("HMM states") + 
                       cowplot::theme_cowplot(font_size = FONT_SIZE) +
                       theme(plot.title = element_text(hjust = 0.5)) +
                       facet_wrap(~state_recode, nrow = N_ROWS) +
                       ggtitle(NULL),
                     x = 0, y = 0,
                     width = 0.333,
                     height = 0.5) +
  cowplot::draw_plot(of_box_sge +
                       cowplot::theme_cowplot(font_size = FONT_SIZE) +
                       theme(plot.title = element_text(hjust = 0.5)) +
                       facet_wrap(~state_recode, nrow = N_ROWS) +
                       ggtitle(NULL),
                     x = 0.333, y = 0,
                     width = 0.333,
                     height = 0.5) +
  cowplot::draw_plot(no_box_sge +
                       cowplot::theme_cowplot(font_size = FONT_SIZE) +
                       theme(plot.title = element_text(hjust = 0.5)) +
                       facet_wrap(~state_recode, nrow = N_ROWS) +
                       ggtitle(NULL),
                     x = 0.666, y = 0,
                     width = 0.333,
                     height = 0.5) +
  cowplot::draw_plot_label(label = c("D", "E", "F"),
                           x = c(0, 0.333, 0.666),
                           y = c(0.5, 0.5, 0.5))

ggsave(POLAR_BOX_DGE_SGE,
       final_dge_sge,
       device = "png",
       width = 16.8,
       height = HEIGHT*2,
       units = "in",
       dpi = 400)

#######################
# Medarkov matrices
#######################

SEC_INT = 2

# DGE

dge_tile_df = dge_df %>% 
  # remove iCab ref fishes (because DGE compares test fishes)
  dplyr::filter(!(line == "iCab" & fish == "ref")) %>% 
  #dplyr::slice_sample(n = 1e6) %>% 
  # add `indiv` column
  tidyr::unite(date, time, quadrant, fish,
               col = "indiv",
               sep = "_") %>%
  # rename and reorder assay
  dplyr::mutate(assay = stringr::str_replace(assay, "_", " "),
                assay = factor(assay, levels = c("open field", "novel object"))) %>% 
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

dge_tile_fig = dge_tile_df %>% 
  ggplot() +
  geom_tile(aes(seconds, indiv, fill = state_recode)) + 
  facet_grid(rows = vars(line), cols = vars(assay), scales = "free") +
  scale_fill_viridis_c() +
  scale_y_discrete(limits = rev) +
  guides(fill = "none") +
  ylab("individual fish") +
  #cowplot::theme_cowplot() +
  #theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

ggsave(TILE_DGE,
       dge_tile_fig,
       device = "png",
       width = 12,
       height = 18,
       units = "in",
       dpi = 400)

# SGE

sge_tile_df = sge_df %>% 
  #dplyr::slice_sample(n = 1e5) %>% 
  # remove iCab test fishes (because SGE compares ref fish)
  dplyr::filter(!(line == "iCab" & fish == "test")) %>% 
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

sge_tile_fig = sge_tile_df %>% 
  ggplot() +
  geom_tile(aes(seconds, indiv, fill = state_recode)) + 
  facet_grid(rows = vars(test_fish), cols = vars(assay), scales = "free") +
  scale_fill_viridis_c(option = "magma") +
  scale_y_discrete(limits = rev) +
  guides(fill = "none") +
  ylab("individual fish") +
  #cowplot::theme_cowplot() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggsave(TILE_SGE,
       sge_tile_fig,
       device = "png",
       width = 12,
       height = 18,
       units = "in",
       dpi = 400)

# Together

tile_dge_sge = cowplot::ggdraw() +
  cowplot::draw_plot(dge_tile_fig +
                       cowplot::theme_cowplot() +
                       theme(strip.background.x = element_blank(),
                             strip.text.y = element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks.y=element_blank(),
                             axis.line.y = element_blank()),
                     x = 0, y = 0,
                     width = 0.5, height = 1) +
  cowplot::draw_plot(sge_tile_fig +
                       cowplot::theme_cowplot() +
                       theme(strip.background = element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks.y=element_blank(),
                             axis.line.y = element_blank()) +
                       ylab(NULL),
                     x = 0.5, y = 0,
                     width = 0.5, height = 1) + 
  cowplot::draw_plot_label(c("A", "B"),
                           x = c(0, 0.5),
                           y = c(1, 1))

ggsave(TILE_DGE_SGE,
       tile_dge_sge,
       device = "png",
       width = 18,
       height = 22,
       units = "in",
       dpi = 400)
