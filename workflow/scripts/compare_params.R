# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)
#library(fishualize)

# Get variables

## Debug
#CONC = list("/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/0.05/dist_angle/10.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/0.05/dist_angle/15.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/0.05/dist_angle/20.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/0.05/dist_angle/5.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/0.1/dist_angle/10.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/0.1/dist_angle/15.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/0.1/dist_angle/20.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/0.1/dist_angle/5.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/0.2/dist_angle/10.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/0.2/dist_angle/15.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/0.2/dist_angle/20.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/0.2/dist_angle/5.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/0.5/dist_angle/10.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/0.5/dist_angle/15.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/0.5/dist_angle/20.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/0.5/dist_angle/5.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/1/dist_angle/10.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/1/dist_angle/15.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/1/dist_angle/20.rds",
#            "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/1/dist_angle/5.rds")
#KW = "/hps/nobackup/birney/users/ian/pilot/kruskal_wallis/out.rds"

## True
CONC = snakemake@input[["conc"]]
KW = snakemake@input[["kw"]]
OUT_PNG = snakemake@output[["png"]]
OUT_PDF = snakemake@output[["pdf"]]

# Read and process

kw = readRDS(KW) %>% 
  dplyr::mutate(INTERVAL = as.numeric(INTERVAL)) %>% 
  dplyr::mutate(N_STATES = as.numeric(N_STATES))

conc = purrr::map_dfr(CONC, readRDS) %>% 
  dplyr::group_by(INTERVAL, N_STATES) %>% 
  # get mean across groups A and B
  dplyr::summarise(MEAN_CONC = mean(CONCORDANCE))

df = dplyr::left_join(kw, conc, by = c("INTERVAL", "N_STATES"))

# Plot

## Get palette
#devtools::source_gist("b68d094672bdde2f0144b755b4721fe6")
#pal = grDevices::colorRampPalette(electro_angler)
pal = colorspace::sequential_hcl(length(unique(df$N_STATES)), palette = "ag_Sunset")

out_plot = df %>% 
  dplyr::mutate(N_STATES = factor(N_STATES, levels = sort(unique(N_STATES)))) %>% 
  dplyr::mutate(INTERVAL = factor(INTERVAL, levels = sort(unique(INTERVAL)))) %>% 
  ggplot() +
    geom_point(aes(MEAN_CONC, SUM_KW_STAT, size = INTERVAL, colour = N_STATES)) +
    theme_bw() +
    guides(size = guide_legend(title = "Interval\n(seconds)"),
           colour = guide_legend(title = "N states")) +
    scale_colour_manual(values = pal) +
    xlab("mean concordance between cross-validated HMM states") +
    ylab("summed Kruskal-Wallis statistic comparing frequency\nof time spent in each HMM state across medaka lines")

ggsave(OUT_PNG,
       out_plot,
       device = "png",
       width = 10,
       height = 8.5,
       units = "in",
       dpi = 400)

ggsave(OUT_PDF,
       out_plot,
       device = "pdf",
       width = 10,
       height = 8.5,
       units = "in",
       dpi = 400)
