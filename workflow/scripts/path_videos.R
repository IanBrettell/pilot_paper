# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)
library(gganimate)
library(wesanderson)

# Get variables

## Debug

IN = "/hps/nobackup/birney/users/ian/pilot/hmm_out/0.08/dist_angle/15.csv"
ASSAY = "open_field"
SAMPLE = "20190613_1307_icab_hdr_R"
INTERVAL = 0.08
DIMS = here::here("config/split_video_dims.csv")

## True

IN = snakemake@input[["hmm"]]
DIMS = snakemake@input[["dims"]]
ASSAY = snakemake@params[["assay"]]
SAMPLE = snakemake@params[["sample"]]
INTERVAL = snakemake@params[["interval"]] %>% 
  as.numeric()
OUT = snakemake@output[[1]]

#######################
# Set plotting parameters
#######################

# Get ref and test

split_samp = stringr::str_split(SAMPLE, pattern = "_", simplify = T)

ref = split_samp[3]
test = split_samp[4]

if (ref == test) {
  line_vec = c("iCab\nref", "iCab\ntest")
  names(line_vec) = c("ref", "test")
  pal = c("#F1BB7B", darker("#F1BB7B", amount = 70))
  names(pal) = line_vec
} else {
  line_vec = c("iCab", "HdrR", "HNI", "Kaga", "HO5")
  names(line_vec) = c("icab", "hdr", "hni", "kaga", "ho5")
  new_pal = grDevices::colorRampPalette(wesanderson::wes_palette("GrandBudapest1"))
  pal = new_pal(5)
  names(pal) = line_vec
  # filter `line_vec` and `pal` for only the lines in the sample
  target_ind = which(names(line_vec) %in% ref  | names(line_vec) %in% test)
  line_vec = line_vec[target_ind]
  pal = pal[target_ind]
}

FPS = 1/INTERVAL

#######################
# Read in data
#######################

# Read in file

df = readr::read_csv(IN) %>% 
  # recode angle to sit between 0 and 360
  dplyr::mutate(angle_recode = ifelse(angle < 0,
                                      180 + (180 + angle),
                                      angle)) %>% 
  # create `sample` column
  tidyr::unite(col = "sample",
               date, time, ref_fish, test_fish, tank_side,
               sep = "_",
               remove = F) %>% 
  # filter for target sample and ref/test
  dplyr::filter(sample == SAMPLE & assay == ASSAY)

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
  dplyr::mutate(state_recode = dplyr::recode(state, !!!recode_vec),
                state_recode = factor(state_recode, levels = recode_vec))

# Recode and order line depending on iCab controls or not

if (ref == test){
  df = df %>% 
    dplyr::mutate(line = dplyr::recode(fish, !!!line_vec)) %>% 
    dplyr::mutate(line = factor(line, levels = line_vec))
} else {
  df = df %>% 
    dplyr::mutate(line = dplyr::case_when(fish == "ref" ~ ref_fish,
                                          fish == "test" ~ test_fish)) %>% 
    # recode
    dplyr::mutate(line = dplyr::recode(line, !!!line_vec)) %>% 
    # factorise to order
    dplyr::mutate(line = factor(line, levels = line_vec))
}

# Read in dims

dims = readr::read_csv(DIMS) %>% 
  dplyr::filter(sample == SAMPLE & assay == ASSAY)

#######################
# Plot
#######################

# Set order of quadrants for facets

QUADRANT_ORDER = c("q2", "q1", "q3", "q4")

# Get max width and height 

wid = dims %>% 
  dplyr::pull(wid) %>% 
  max()
hei = dims %>% 
  dplyr::pull(hei) %>% 
  max()

# Get total height and width in pixels (to match with the raw videos)

TOT_WID = dims %>%
  dplyr::filter(quadrant %in% c("q1", "q2")) %>% 
  dplyr::pull(wid) %>% 
  sum()

TOT_HEI = dims %>%
  dplyr::filter(quadrant %in% c("q1", "q4")) %>% 
  dplyr::pull(hei) %>% 
  sum()

# Get total number of frames

N_FRAMES = df %>% 
  dplyr::distinct(seconds) %>% 
  nrow()

out = df %>% 
  #dplyr::filter(seconds < 30) %>% 
  dplyr::mutate(quadrant = factor(quadrant, levels = QUADRANT_ORDER)) %>% 
  ggplot() +
  geom_path(aes(x, y, colour = line)) +
  facet_wrap(~quadrant, nrow = 2) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(limits = c(0,wid)) +
  scale_y_reverse(limits = c(hei,0)) +
  guides(colour = "none") +
  cowplot::theme_cowplot(font_size = 10) +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()
        ) +
  gganimate::transition_reveal(seconds, keep_last = T)

gganimate::animate(
  out,
  width = 9,
  height = 9,
  units = "in",
  res = 400,
  nframes = N_FRAMES,
  fps = FPS,
  renderer = gganimate::av_renderer(OUT),
  device = "png"
)
