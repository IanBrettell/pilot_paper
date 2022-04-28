# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(ggbeeswarm)
library(devtools)
library(cowplot)
library(ggpubr)

## Debug
#IN = "/hps/nobackup/birney/users/ian/pilot/hmm_out/0.05/dist_angle/18.csv"

## True
IN = snakemake@input[[1]]
OUT = snakemake@output[["fig"]]

# Get lighter/darker functions

devtools::source_gist("c5015ee666cdf8d9f7e25fa3c8063c99")

# Read in file

df = readr::read_csv(IN) 

df_control = df %>% 
  # Filter for only iCabs paired with iCabs
  dplyr::filter(test_fish == "icab") %>% 
  # Get individual
  tidyr::unite(date, time, quadrant, fish,
               col = "indiv", 
               remove = F) %>% 
  # Group by assay and individual to get mean speed
  dplyr::group_by(assay, date, time, quadrant, indiv) %>% 
  # Calculate mean speed
  dplyr::summarise(mean_speed = mean(distance)) %>% 
  dplyr::ungroup() %>% 
  # Make date a factor
  dplyr::mutate(date = as.factor(date))

###################
# Linear model
###################

model_fit = lm(mean_speed ~ assay + date + time + quadrant, data = df_control)
summary(model_fit)

###################
# Kruskal-Wallis tests
###################

#kw = tibble::tibble("covariate" = c("date", "time", "quadrant")) %>% 
#  dplyr::mutate("p_value" = kruskal.test(df_control$mean_speed, g = df_control[[covariate]])$p.value)
#  
#
#
#kw = tibble::tibble("covariate" = c("date", "time", "quadrant"),
#                    "p" = c(kruskal.test(df_control$mean_speed, g = df_control$date)$p.value,
#                             kruskal.test(df_control$mean_speed, g = df_control$time)$p.value,
#                             kruskal.test(df_control$mean_speed, g = df_control$quadrant)$p.value)) %>% 
#  rstatix::adjust_pvalue(method = "fdr") %>% 
#  dplyr::mutate(p.adj = signif(p.adj, digits = 3)) %>% 
#  rstatix::add_significance(p.col = "p.adj") %>% 
#  # paste p-value and significance together
#  dplyr::mutate(p_final = dplyr::case_when(p.adj.signif == "ns" ~ paste("p =", p.adj),
#                                           TRUE ~ paste("p =", p.adj, p.adj.signif)))

# Plots

## Date

date_pal = colorspace::sequential_hcl(length(unique(df_control$date)),
                                      palette = "OrYel")
  
#kw_date = kruskal.test(df_control$mean_speed, g = df_control$date)$p.value %>% 
#  signif(., digits = 3)

date_fig = df_control %>% 
  dplyr::mutate(date = factor(date, levels = sort(unique(date)))) %>% 
  ggplot(aes(date, mean_speed, fill = date)) +
  geom_violin(aes(colour = date)) +
  geom_boxplot(aes(colour = date),
               width = 0.3) +
  ggbeeswarm::geom_beeswarm(colour = "#3B1F2B", alpha = 0.8) +
  #geom_text(x = "20190611", y = 3.15, label = paste("p =", kw_date)) +
  scale_fill_manual(values = date_pal) +
  scale_colour_manual(values = darker(date_pal, amount = 100)) +
  cowplot::theme_cowplot() +
  guides(fill = "none",
         colour = "none") +
  ylab("mean speed") +
  ggpubr::stat_compare_means(label.y = 3.15)

## Time

time_pal = colorspace::sequential_hcl(length(unique(df_control$time)),
                                      palette = "ag_GrnYl")

#kw_time = kruskal.test(df_control$mean_speed, g = df_control$time)$p.value %>% 
#  signif(., digits = 3)

time_fig = df_control %>% 
  dplyr::mutate(time = factor(time)) %>% 
  ggplot(aes(time, mean_speed, fill = time)) +
  geom_violin(aes(colour = time)) +
  geom_boxplot(aes(colour = time),
               width = 0.2) +
  ggbeeswarm::geom_beeswarm(colour = "#3B1F2B", alpha = 0.8) +
  #geom_text(x = "I", y = 3.15, label = paste("p =", kw_time)) +
  scale_fill_manual(values = time_pal) +
  scale_colour_manual(values = darker(time_pal, amount = 100)) +
  cowplot::theme_cowplot() +
  guides(colour = "none",
         fill = "none") +
  ylab("mean speed") +
  ggpubr::stat_compare_means(label.x = "952", label.y = 3.15)

## Quadrant
  
#quad_pal = colorspace::sequential_hcl(length(unique(df_control$quadrant)),
#                                      palette = "magenta")
quad_pal = scales::hue_pal()(4)

#kw_quad = kruskal.test(df_control$mean_speed, g = df_control$quadrant)$p.value %>% 
#  signif(., digits = 3)

quad_fig = df_control %>% 
  dplyr::mutate(quadrant = dplyr::recode(quadrant,
                                         "q1" = "I",
                                         "q2" = "II",
                                         "q3" = "III",
                                         "q4" = "IV")) %>% 
  dplyr::mutate(quadrant = factor(quadrant, levels = c("I", "II", "III", "IV"))) %>% 
  ggplot(aes(quadrant, mean_speed, fill = quadrant)) +
  geom_violin(aes(colour = quadrant)) +
  geom_boxplot(aes(colour = quadrant),
               width = 0.3) +
  ggbeeswarm::geom_beeswarm(colour = "#3B1F2B", alpha = 0.8) +
  #geom_text(x = "I", y = 3.15, label = paste("p =", kw_quad)) +
  scale_fill_manual(values = quad_pal) +
  scale_colour_manual(values = darker(quad_pal, amount = 100)) +
  cowplot::theme_cowplot() +
  guides(colour = "none",
         fill = "none") +
  ylab("mean speed") +
  ggpubr::stat_compare_means(label.y = 3.15)

###################
# Compile into final
###################

final_fig = cowplot::ggdraw() +
  cowplot::draw_plot(date_fig,
                     x = 0, y = 0.5,
                     width = 0.55, height = 0.5) +
  cowplot::draw_plot(quad_fig,
                     x = 0.55, y = 0.5,
                     width = 0.45, height = 0.5) +
  cowplot::draw_plot(time_fig,
                     x = 0, y = 0,
                     width = 1, height = 0.5) +
  cowplot::draw_plot_label(c("A", "B", "C"),
                           x = c(0, 0.55, 0),
                           y = c(1, 1, 0.5))

ggsave(OUT,
       final_fig,
       device = "png",
       width = 12,
       height = 12,
       units = "in",
       dpi = 400)
