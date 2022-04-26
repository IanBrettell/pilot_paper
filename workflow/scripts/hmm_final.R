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
IN = "/hps/nobackup/birney/users/ian/pilot/hmm_out/0.08/dist_angle/12.csv"
N_STATES = 12
VARIABLES = "distance and angle of travel"
INTERVAL = 0.08

#######################
# Read in data
#######################

# Get number of rows (for plotting) based on number of states

if (N_STATES %% 5 == 0){
  N_ROWS = N_STATES / 5
} else {
  N_ROWS = 3
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

#######################
# Polar plot
#######################

# Set title
TITLE = paste("N states: ",
              N_STATES,
              "\nVariables: ",
              VARIABLES,
              "\nInterval: ",
              INTERVAL, " seconds",
              sep = "")

out_plot = df %>% 
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

dge_df %>% 
  # get frequency of time spent in each state for each line
  dplyr::group_by(line) %>% 
  dplyr::count(line, state_recode) %>% 
  dplyr::add_count(line, wt = n, name = "nn") %>% 
  dplyr::mutate(state_freq = n / nn) %>% 
  ggplot() +
  geom_col(aes(state_recode, state_freq, fill = state_recode)) +
  facet_grid(rows = vars(line)) +
  scale_fill_viridis_c() +
  theme_bw() +
  scale_x_continuous(breaks = unique(dge_df$state_recode)) +
  guides(fill = "none") +
  xlab("HMM state") +
  ylab("Proportion of time spent in HMM state")


#######################
# Box plots
#######################  

# Get proportions of time spent in each state

state_freq_df = dge_df %>% 
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

kw_dge_2 = state_freq_df %>% 
  dplyr::group_by(assay, state_recode) %>% 
  rstatix::kruskal_test(state_freq ~ line) %>% 
  rstatix::adjust_pvalue(method = "fdr") %>% 
  dplyr::mutate(p.adj = signif(p.adj, digits = 3)) %>% 
  rstatix::add_significance(p.col = "p.adj") %>% 
  # paste p-value and significance together
  dplyr::mutate(p_final = paste(italic(p), "=", p.adj, p.adj.signif))


# Plot boxplots

## Open field

ASSAY = "open_field"
state_freq_df %>% 
  dplyr::filter(assay == ASSAY) %>% 
  ggplot() +
    geom_boxplot(aes(line, state_freq, fill = line), notch = T) +
    #ggbeeswarm::geom_beeswarm(aes(line, state_freq), size = 0.5) +
    #ggpubr::stat_pvalue_manual(data = kw_dge_2, label = "p.adj", remove.bracket = T, y.position = 0.8) +
    geom_text(data = kw_dge_2 %>% 
                dplyr::filter(assay == ASSAY),
              aes(x = "iCab", y = 0.8, label = p.adj),
              size = 3) +
    scale_fill_manual(values = pal) +
    theme_bw() +
    facet_wrap(~state_recode, nrow = 3) +
    guides(fill = "none") +
    ylab("HMM state frequency") +
    ggtitle(stringr::str_replace(ASSAY, "_", " ")) +
    theme(plot.title = element_text(hjust = 0.5))

## Novel object

state_freq_df %>% 
  dplyr::filter(assay == "novel_object") %>% 
  ggplot(aes(line, state_freq)) +
  geom_boxplot(aes(fill = line), notch = T) +
  #ggbeeswarm::geom_beeswarm(aes(line, state_freq), size = 0.5) +
  ggpubr::stat_compare_means(method = "kruskal.test") +
  scale_fill_manual(values = pal) +
  theme_bw() +
  facet_wrap(~state_recode, nrow = N_ROWS)
