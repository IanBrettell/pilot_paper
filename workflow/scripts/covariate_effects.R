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
IN = "/hps/nobackup/birney/users/ian/pilot/hmm_out/0.05/dist_angle/18.csv"

## True
IN = snakemake@input[[1]]
OUT_OF = snakemake@output[["of"]]
OUT_NO = snakemake@output[["no"]]

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
  dplyr::group_by(assay, date, time, quadrant, tank_side, indiv) %>% 
  # Calculate mean speed
  dplyr::summarise(mean_speed = mean(distance)) %>% 
  dplyr::ungroup() %>% 
  # Make date a factor
  dplyr::mutate(date = as.factor(date))

###################
# Linear model
###################

model_fit = lm(mean_speed ~ assay + date + time + quadrant + tank_side, data = df_control)
summary(model_fit)

###################
# ANOVA model
###################

aov_fit = aov(mean_speed ~ assay + date + time + quadrant + tank_side, data = df_control)
summary(aov_fit)

# Without day 1

aov_fit_noday1 = aov(mean_speed ~ assay + date + time + quadrant + tank_side,
                     data = df_control %>% 
                       dplyr::filter(date != "20190611"))
summary(aov_fit_noday1)

# On all simultaneously

aov_df = df_control %>% 
  dplyr::group_by(assay) %>% 
  tidyr::nest() %>%
  dplyr::mutate(model = purrr::map(data, ~aov(
    mean_speed ~ date + time + quadrant + tank_side,
    data = .))) %>%
  dplyr::select(-data) %>% 
  dplyr::mutate(model_tidy = purrr::map(model, broom::tidy)) %>%
  tidyr::unnest(model_tidy) %>% 
  rstatix::add_significance(p.col = "p.value") %>% 
  # remove model
  dplyr::select(-model) %>% 
  # reduce to 3 digits
  dplyr::mutate(p.value = signif(p.value, digits = 3)) %>% 
  # paste p-value with significance
  dplyr::mutate(p_final = dplyr::case_when(p.value.signif == "ns" ~ paste("p =", p.value),
                                           TRUE ~ paste("p =", p.value, p.value.signif)))

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


#########################
# Plot OF
#########################


df_of = df_control %>% 
  dplyr::filter(assay == "open_field")

## Date

date_pal = colorspace::sequential_hcl(length(unique(df_control$date)),
                                      palette = "OrYel")
  
#kw_date = kruskal.test(df_control$mean_speed, g = df_control$date)$p.value %>% 
#  signif(., digits = 3)

date_fig_of = df_of %>% 
  dplyr::mutate(date = factor(date, levels = sort(unique(date)))) %>% 
  #ggplot(aes(date, mean_speed, fill = date)) +
  ggplot() +
  geom_violin(aes(date, mean_speed, fill = date, colour = date)) +
  geom_boxplot(aes(date, mean_speed, fill = date, colour = date),
               width = 0.25) +
  ggbeeswarm::geom_beeswarm(aes(date, mean_speed),
                            colour = "#3B1F2B", alpha = 0.8) +
  geom_text(data = aov_df %>% 
              dplyr::filter(assay == "open_field" & term == "date"),
            aes(x = "20190611", y = 3, label = p_final)) +
  #geom_text(x = "20190611", y = 3.15, label = paste("p =", kw_date)) +
  scale_fill_manual(values = date_pal) +
  scale_colour_manual(values = darker(date_pal, amount = 100)) +
  cowplot::theme_cowplot() +
  guides(fill = "none",
         colour = "none") +
  ylab("mean speed") 

## Time

time_pal = colorspace::sequential_hcl(length(unique(df_control$time)),
                                      palette = "ag_GrnYl")

#kw_time = kruskal.test(df_control$mean_speed, g = df_control$time)$p.value %>% 
#  signif(., digits = 3)

time_fig_of = df_of %>% 
  dplyr::mutate(time = factor(time)) %>% 
  ggplot() +
  geom_violin(aes(time, mean_speed, fill = time, colour = time)) +
  geom_boxplot(aes(time, mean_speed, fill = time, colour = time),
               width = 0.15) +
  ggbeeswarm::geom_beeswarm(aes(time, mean_speed),
                            colour = "#3B1F2B", alpha = 0.8) +
  geom_text(data = aov_df %>% 
              dplyr::filter(assay == "open_field" & term == "time"),
            aes(x = "952", y = 3, label = p_final)) +
  scale_fill_manual(values = time_pal) +
  scale_colour_manual(values = darker(time_pal, amount = 100)) +
  cowplot::theme_cowplot() +
  guides(colour = "none",
         fill = "none") +
  ylab("mean speed") 

## Quadrant
  
#quad_pal = colorspace::sequential_hcl(length(unique(df_control$quadrant)),
#                                      palette = "magenta")
quad_pal = scales::hue_pal()(4)

#kw_quad = kruskal.test(df_control$mean_speed, g = df_control$quadrant)$p.value %>% 
#  signif(., digits = 3)

quad_fig_of = df_of %>% 
  dplyr::mutate(quadrant = dplyr::recode(quadrant,
                                         "q1" = "I",
                                         "q2" = "II",
                                         "q3" = "III",
                                         "q4" = "IV")) %>% 
  dplyr::mutate(quadrant = factor(quadrant, levels = c("I", "II", "III", "IV"))) %>% 
  #ggplot(aes(quadrant, mean_speed, fill = quadrant)) +
  ggplot() +
  geom_violin(aes(quadrant, mean_speed, fill = quadrant, colour = quadrant)) +
  geom_boxplot(aes(quadrant, mean_speed, fill = quadrant, colour = quadrant),
               width = 0.3) +
  ggbeeswarm::geom_beeswarm(aes(quadrant, mean_speed),
                            colour = "#3B1F2B", alpha = 0.8) +
  geom_text(data = aov_df %>% 
              dplyr::filter(assay == "open_field" & term == "quadrant"),
            aes(x = "I", y = 3.15, label = p_final)) +
  scale_fill_manual(values = quad_pal) +
  scale_colour_manual(values = darker(quad_pal, amount = 100)) +
  cowplot::theme_cowplot() +
  guides(colour = "none",
         fill = "none") +
  ylab("mean speed") 

## Tank side

tank_pal = colorspace::diverging_hcl(length(unique(df_control$tank_side)),
                                      palette = "Red-Green")

tank_fig_of = df_of %>% 
  dplyr::mutate(tank_side = factor(tank_side, levels = c("L", "R"))) %>% 
  #ggplot(aes(tank_side, mean_speed, fill = tank_side)) +]
  ggplot() +
  geom_violin(aes(tank_side, mean_speed, fill = tank_side, colour = tank_side)) +
  geom_boxplot(aes(tank_side, mean_speed, fill = tank_side, colour = tank_side),
               width = 0.25) +
  ggbeeswarm::geom_beeswarm(aes(tank_side, mean_speed),
                            colour = "#3B1F2B", alpha = 0.8) +  
  geom_text(data = aov_df %>% 
              dplyr::filter(assay == "open_field" & term == "tank_side"),
            aes(x = "L", y = 3.2, label = p_final)) +
  scale_fill_manual(values = tank_pal) +
  scale_colour_manual(values = darker(tank_pal, amount = 100)) +
  cowplot::theme_cowplot() +
  guides(colour = "none",
         fill = "none") +
  ylab("mean speed") +
  xlab("tank side")


## Compile

final_of = cowplot::ggdraw() +
  cowplot::draw_plot(date_fig_of,
                     x = 0, y = 0.5,
                     width = 0.55, height = 0.5) + 
  cowplot::draw_plot(quad_fig_of,
                     x = 0.55, y = 0.5,
                     width = 0.45, height = 0.5) + 
  cowplot::draw_plot(time_fig_of,
                     x = 0, y = 0,
                     width = 0.6, height = 0.5) + 
  cowplot::draw_plot(tank_fig_of,
                     x = 0.6, y = 0,
                     width = 0.4, height = 0.5) +
  cowplot::draw_plot_label(c("A", "B", "C", "D"),
                           x = c(0, 0.55, 0, 0.6),
                           y = c(1, 1, 0.5, 0.5))
  

ggsave(OUT_OF,
       final_of,
       device = "png",
       width = 15,
       height = 9,
       units = "in",
       dpi = 400)


#########################
# Plot NO
#########################

ASSAY = "novel_object"

df_no = df_control %>% 
  dplyr::filter(assay == ASSAY)

## Date

date_pal = colorspace::sequential_hcl(length(unique(df_control$date)),
                                      palette = "OrYel")

#kw_date = kruskal.test(df_control$mean_speed, g = df_control$date)$p.value %>% 
#  signif(., digits = 3)

date_fig_no = df_no %>% 
  dplyr::mutate(date = factor(date, levels = sort(unique(date)))) %>% 
  #ggplot(aes(date, mean_speed, fill = date)) +
  ggplot() +
  geom_violin(aes(date, mean_speed, fill = date, colour = date)) +
  geom_boxplot(aes(date, mean_speed, fill = date, colour = date),
               width = 0.25) +
  ggbeeswarm::geom_beeswarm(aes(date, mean_speed),
                            colour = "#3B1F2B", alpha = 0.8) +
  geom_text(data = aov_df %>% 
              dplyr::filter(assay == ASSAY & term == "date"),
            aes(x = "20190611", y = 3, label = p_final)) +
  #geom_text(x = "20190611", y = 3.15, label = paste("p =", kw_date)) +
  scale_fill_manual(values = date_pal) +
  scale_colour_manual(values = darker(date_pal, amount = 100)) +
  cowplot::theme_cowplot() +
  guides(fill = "none",
         colour = "none") +
  ylab("mean speed") 

## Time

time_pal = colorspace::sequential_hcl(length(unique(df_control$time)),
                                      palette = "ag_GrnYl")

#kw_time = kruskal.test(df_control$mean_speed, g = df_control$time)$p.value %>% 
#  signif(., digits = 3)

time_fig_no = df_no %>% 
  dplyr::mutate(time = factor(time)) %>% 
  ggplot() +
  geom_violin(aes(time, mean_speed, fill = time, colour = time)) +
  geom_boxplot(aes(time, mean_speed, fill = time, colour = time),
               width = 0.15) +
  ggbeeswarm::geom_beeswarm(aes(time, mean_speed),
                            colour = "#3B1F2B", alpha = 0.8) +
  geom_text(data = aov_df %>% 
              dplyr::filter(assay == ASSAY & term == "time"),
            aes(x = "952", y = 3, label = p_final)) +
  scale_fill_manual(values = time_pal) +
  scale_colour_manual(values = darker(time_pal, amount = 100)) +
  cowplot::theme_cowplot() +
  guides(colour = "none",
         fill = "none") +
  ylab("mean speed") 

## Quadrant

#quad_pal = colorspace::sequential_hcl(length(unique(df_control$quadrant)),
#                                      palette = "magenta")
quad_pal = scales::hue_pal()(4)

#kw_quad = kruskal.test(df_control$mean_speed, g = df_control$quadrant)$p.value %>% 
#  signif(., digits = 3)

quad_fig_no = df_no %>% 
  dplyr::mutate(quadrant = dplyr::recode(quadrant,
                                         "q1" = "I",
                                         "q2" = "II",
                                         "q3" = "III",
                                         "q4" = "IV")) %>% 
  dplyr::mutate(quadrant = factor(quadrant, levels = c("I", "II", "III", "IV"))) %>% 
  #ggplot(aes(quadrant, mean_speed, fill = quadrant)) +
  ggplot() +
  geom_violin(aes(quadrant, mean_speed, fill = quadrant, colour = quadrant)) +
  geom_boxplot(aes(quadrant, mean_speed, fill = quadrant, colour = quadrant),
               width = 0.3) +
  ggbeeswarm::geom_beeswarm(aes(quadrant, mean_speed),
                            colour = "#3B1F2B", alpha = 0.8) +
  geom_text(data = aov_df %>% 
              dplyr::filter(assay == ASSAY & term == "quadrant"),
            aes(x = "I", y = 3.15, label = p_final)) +
  scale_fill_manual(values = quad_pal) +
  scale_colour_manual(values = darker(quad_pal, amount = 100)) +
  cowplot::theme_cowplot() +
  guides(colour = "none",
         fill = "none") +
  ylab("mean speed") 

## Tank side

tank_pal = colorspace::diverging_hcl(length(unique(df_control$tank_side)),
                                     palette = "Red-Green")

tank_fig_no = df_no %>% 
  dplyr::mutate(tank_side = factor(tank_side, levels = c("L", "R"))) %>% 
  #ggplot(aes(tank_side, mean_speed, fill = tank_side)) +]
  ggplot() +
  geom_violin(aes(tank_side, mean_speed, fill = tank_side, colour = tank_side)) +
  geom_boxplot(aes(tank_side, mean_speed, fill = tank_side, colour = tank_side),
               width = 0.25) +
  ggbeeswarm::geom_beeswarm(aes(tank_side, mean_speed),
                            colour = "#3B1F2B", alpha = 0.8) +  
  geom_text(data = aov_df %>% 
              dplyr::filter(assay == ASSAY & term == "tank_side"),
            aes(x = "L", y = 3.2, label = p_final)) +
  scale_fill_manual(values = tank_pal) +
  scale_colour_manual(values = darker(tank_pal, amount = 100)) +
  cowplot::theme_cowplot() +
  guides(colour = "none",
         fill = "none") +
  ylab("mean speed") +
  xlab("tank side")


## Compile

final_no = cowplot::ggdraw() +
  cowplot::draw_plot(date_fig_no,
                     x = 0, y = 0.5,
                     width = 0.55, height = 0.5) + 
  cowplot::draw_plot(quad_fig_no,
                     x = 0.55, y = 0.5,
                     width = 0.45, height = 0.5) + 
  cowplot::draw_plot(time_fig_no,
                     x = 0, y = 0,
                     width = 0.6, height = 0.5) + 
  cowplot::draw_plot(tank_fig_no,
                     x = 0.6, y = 0,
                     width = 0.4, height = 0.5) +
  cowplot::draw_plot_label(c("A", "B", "C", "D"),
                           x = c(0, 0.55, 0, 0.6),
                           y = c(1, 1, 0.5, 0.5))


ggsave(OUT_NO,
       final_no,
       device = "png",
       width = 15,
       height = 9,
       units = "in",
       dpi = 400)



