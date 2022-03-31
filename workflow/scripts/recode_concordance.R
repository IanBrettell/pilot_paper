# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)

# Get variables

## Debug
IN_FILE = "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_out/0.1/dist_angle/15.csv"
N_STATES = 15

## True
IN_FILE = snakemake@input[[1]]
OUT_CSV = snakemake@output[["csv"]]
OUT_PNG = snakemake@output[["png"]]
OUT_PDF = snakemake@output[["pdf"]]
N_STATES = snakemake@params[["n_states"]] %>% 
  as.numeric()


# Get number of rows (for plotting) based on number of states

N_ROWS = N_STATES / 5

# Read in file

df = readr::read_csv(IN_FILE) %>% 
  # recode angle to sit between 0 and 360
  dplyr::mutate(angle_recode = ifelse(angle < 0,
                                      180 + (180 + angle),
                                      angle))

# Rank states by mean distance

rank_df = df %>% 
  #dplyr::slice_sample(n = 1e4) %>% 
  tidyr::pivot_longer(cols = c(train_self, train_other),
                      names_to = "training_data",
                      names_prefix = "train_",
                      values_to = "state") %>% 
  dplyr::mutate(log10_dist = log10(distance + 1))

# Run PCA on distance and angle
pca_out = prcomp(rank_df %>% 
                   dplyr::select(log10_dist, angle),
                 center = TRUE, scale. = T)
pc1 = pca_out$x[, 1]
  
# Add to rank_df

rank_df = rank_df %>% 
  dplyr::mutate(pc1 = pc1) %>% 
  dplyr::group_by(group, training_data, state) %>% 
  dplyr::summarise(mean_pc1 = mean(pc1)) %>% 
  # order by mean distance
  dplyr::arrange(mean_pc1, .by_group = T) %>% 
  # recode state by ranking
  dplyr::mutate(state_recode = dplyr::row_number()) %>% 
  dplyr::ungroup()

# Bind to original DF

df_new = dplyr::left_join(df %>% 
                            tidyr::pivot_longer(cols = c(train_self, train_other),
                                                names_to = "training_data",
                                                names_prefix = "train_",
                                                values_to = "state") , 
                          rank_df ,
                          by = c("group", "training_data", "state"))

# Plot

## Split by group/training_data combinations
df_split = df_new %>% 
  # order `training_data`
  dplyr::mutate(training_data = factor(training_data, levels = c("self", "other"))) %>% 
  split(., ~ group + training_data)

counter = 0
out_plots = lapply(df_split, function(DATA){
  counter <<- counter + 1
  
  # get plot title 
  title = names(df_split)[counter] %>% 
    stringr::str_split(pattern = "\\.", simplify = T)
  title = paste("Group: ", title[1], "\nTraining data: ", title[2], sep = "")
  # Plot
  out_plot = DATA %>% 
    # select random sample of 1e5 rows
    dplyr::slice_sample(n = 1e5) %>% 
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
    ggtitle(title)
  })

# Compile 4 plots into single plot

final_plot = cowplot::plot_grid(plotlist = out_plots,
                                nrow = 2,
                                ncol = 2,
                                axis = "tblr")

# Save .png
ggsave(OUT_PNG,
       final_plot,
       device = "png",
       width = 19.416,
       height = 12,
       units = "in",
       dpi = 400)

# Save .pdf
ggsave(OUT_PDF,
       final_plot,
       device = "pdf",
       width = 19.416,
       height = 12,
       units = "in",
       dpi = 400)


# Pivot wider and save to .csv

df_new %>% 
  # remove `state` column because it prevents proper widening (as it differs based on training data)
  dplyr::select(-state) %>% 
  tidyr::pivot_wider(names_from = training_data,
                     names_prefix = "trained_",
                     values_from = "state_recode") %>% 
  readr::write_csv(OUT_CSV)
