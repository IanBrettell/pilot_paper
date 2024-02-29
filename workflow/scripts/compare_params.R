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
CONC = as.list(list.files("/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/",full.names = T, recursive = T))
KW = "/hps/nobackup/birney/users/ian/pilot/kruskal_wallis/out.rds"
OUT_PNG = here::here("book/figs/compare_params/compare_params.png")
OUT_PDF = here::here("book/figs/compare_params/compare_params.pdf")

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

plot_df = df %>% 
  dplyr::mutate(N_STATES = factor(N_STATES, levels = sort(unique(N_STATES)))) %>% 
  dplyr::mutate(INTERVAL = factor(INTERVAL, levels = sort(unique(INTERVAL)))) %>% 
  # make new column combining `INTERVAL` and `N_STATES`
  tidyr::unite(INTERVAL, N_STATES,
               col = "INT_STATE",
               sep = ";",
               remove = F)

# get middle of x-axis
mid_x = ((max(plot_df$MEAN_CONC) -min(plot_df$MEAN_CONC)) / 2) + min(plot_df$MEAN_CONC)
x_len = (max(plot_df$MEAN_CONC) -min(plot_df$MEAN_CONC))/5
x_r = mid_x + x_len
x_l = mid_x - x_len

# get middle of y-axis
mid_y = ((max(plot_df$SUM_KW_STAT) -min(plot_df$SUM_KW_STAT)) / 2) + min(plot_df$SUM_KW_STAT)
y_len = (max(plot_df$SUM_KW_STAT) -min(plot_df$SUM_KW_STAT))/5
y_t = mid_y + y_len
y_b = mid_y - y_len

# get ranges for diagonal
len_diag = 0.15
start_x = (((max(plot_df$MEAN_CONC) -min(plot_df$MEAN_CONC)))) * (1 - len_diag) + min(plot_df$MEAN_CONC)
start_y = (((max(plot_df$SUM_KW_STAT) -min(plot_df$SUM_KW_STAT)))) * (1 - len_diag) + min(plot_df$SUM_KW_STAT)

out_plot = plot_df %>% 
  ggplot(aes(MEAN_CONC, SUM_KW_STAT)) +
  geom_point(aes(size = INTERVAL, colour = N_STATES),
             alpha = 0.8) +
  ggrepel::geom_text_repel(aes(label = INT_STATE),
                           size = 2,
  ) +
  cowplot::theme_cowplot() +
  guides(size = guide_legend(title = "Interval\n(seconds)"),
         colour = guide_legend(title = "N states")) +
  scale_colour_manual(values = pal) +
  xlab("mean concordance between cross-validated HMM states") +
  ylab("summed Kruskal-Wallis statistic comparing frequency\nof time spent in each HMM state across medaka lines") +
  # add annotations
  ## arrows
  annotate("segment", x=0.33, y=y_b, xend=0.33, yend=y_t,
           colour = "#520026",
           arrow=arrow(length=unit(0.3, "cm")),
           size = 1.2) +
  annotate("segment", x=x_l, y=-55, xend=x_r, yend=-55,
           colour = "#520026",
           arrow=arrow(length=unit(0.3, "cm")),
           size = 1.2) +
  ## text for arrows
  annotate("text", x=0.33, y=((y_t - y_b)/2) + y_b,
           colour = "#520026",
           label="better at distinguishing\nbetween lines",
           angle = 90,
           size = 4.5) +
  annotate("text", x=((x_r - x_l)/2) + x_l, y = -55,
           colour = "#520026",
           label="less overfitting",
           vjust = -1,
           size = 4.5) +
  # diagonal arrow
  annotate("segment", x = start_x, y = start_y, xend = max(plot_df$MEAN_CONC), yend = max(plot_df$SUM_KW_STAT),
           colour = "#C9005E",
           arrow=arrow(length=unit(0.3, "cm")),
           size = 1.4) +
  coord_cartesian(xlim = c(min(plot_df$MEAN_CONC),
                           max(plot_df$MEAN_CONC)),
                  ylim = c(min(plot_df$SUM_KW_STAT),
                           max(plot_df$SUM_KW_STAT)),
                  clip="off") +
  theme(plot.margin = unit(c(0,0,3.5,5), "lines"))

#out_plot

ggsave(OUT_PNG,
       out_plot,
       device = "png",
       width = 11.2,
       height = 9.5,
       units = "in",
       dpi = 400)

ggsave(OUT_PDF,
       out_plot,
       device = "pdf",
       width = 10,
       height = 8.5,
       units = "in",
       dpi = 400)
