# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)
library(wesanderson)

# Get variables

## Debug
IN = "/hps/nobackup/birney/users/ian/pilot/hmm_out/0.08/dist_angle/15.csv"
DIMS = here::here("config/split_video_dims.csv")
SAMPLE = "20190616_1717_icab_icab_L"
SSHOTS = list("/hps/software/users/birney/ian/repos/pilot_paper/results/split_coord_images/open_field/20190611_1331_icab_icab_R.png",
              "/hps/software/users/birney/ian/repos/pilot_paper/results/split_coord_images/novel_object/20190611_1331_icab_icab_R.png") %>% 
  unlist() %>% 
  # Put into correct order
  sort(.,decreasing = T)

## True
IN = snakemake@input[["data"]]
DIMS = snakemake@input[["dims"]]
SAMPLE = snakemake@params[["sample"]]
OUT = snakemake@output[["paths"]]

###################
# Set parameters
###################

MAX_SECONDS = 180

# Get `darker()` and `lighter()` functions

devtools::source_gist("c5015ee666cdf8d9f7e25fa3c8063c99")

# Get ref and test

split_samp = stringr::str_split(SAMPLE, pattern = "_", simplify = T)

ref = split_samp[3]
test = split_samp[4]

if (ref == test) {
  line_vec = c("iCab ref", "iCab test")
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

###################
# Read in files
###################

dims = readr::read_csv(DIMS)

df = readr::read_csv(IN) %>% 
  # create `sample` variable
  tidyr::unite(col = "sample",
               date, time, ref_fish, test_fish, tank_side,
               sep = "_",
               remove = F) %>% 
  # filter for target sample
  dplyr::filter(sample == SAMPLE)

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

# Arrange by assay, frame

df = df %>% 
  dplyr::arrange(assay, line, frame)

###################
# Plot
###################

ASSAYS = c("open_field", "novel_object"); names(ASSAYS) = ASSAYS
QUADRANTS = c("q1", "q2", "q3", "q4"); names(QUADRANTS) = QUADRANTS

# Loop over assay
plot_list = purrr::map(ASSAYS, function(ASSAY){
  
  # Loop over quadrant
  quad_list = purrr::map(QUADRANTS, function(QUADRANT){
    
    # Get dimensions
    target_row = dims %>% 
      dplyr::filter(sample == SAMPLE & assay == ASSAY & quadrant == QUADRANT) 
    wid = target_row %>% 
      dplyr::pull(wid)
    hei = target_row %>% 
      dplyr::pull(hei)
      
    # Plot
    out = df %>% 
      dplyr::filter(assay == ASSAY & quadrant == QUADRANT & seconds > 0 & seconds <= MAX_SECONDS) %>% 
      ggplot() +
      geom_path(aes(x, y, colour = line)) +
      scale_colour_manual(values = pal) +
      coord_fixed() +
      cowplot::theme_cowplot() +
      scale_y_reverse(limits = c(hei, 0)) +
      xlim(0, wid)
    
    return(out)
  })

})

# Compose OF
pgrid_of = cowplot::plot_grid(plot_list$open_field$q2 +
                                theme(legend.position="none",
                                      axis.title = element_blank()),
                              plot_list$open_field$q1 +
                                theme(legend.position="none",
                                      axis.title = element_blank()),
                              plot_list$open_field$q3 +
                                theme(legend.position="none",
                                      axis.title = element_blank()),
                              plot_list$open_field$q4 +
                                theme(legend.position="none",
                                      axis.title = element_blank())
                              )

# get legend
legend = cowplot::get_legend(plot_list$open_field$q1)
#y_axis = cowplot::get_y_axis(plot_list$open_field$q1)
x_lab = cowplot::get_plot_component(plot_list$open_field$q1, "xlab-b")
y_lab = cowplot::get_plot_component(plot_list$open_field$q1, "ylab-l")

pgrid_no = cowplot::plot_grid(plot_list$novel_object$q2 +
                                theme(legend.position="none",
                                      axis.title = element_blank()),
                              plot_list$novel_object$q1 +
                                theme(legend.position="none",
                                      axis.title = element_blank()),
                              plot_list$novel_object$q3 +
                                theme(legend.position="none",
                                      axis.title = element_blank()),
                              plot_list$novel_object$q4+
                                theme(legend.position="none",
                                      axis.title = element_blank())
                              )

final = cowplot::plot_grid(y_lab, pgrid_of, pgrid_no, legend, ncol = 4, rel_widths = c(0.1, 1, 1, 0.1))

plot_grid(final, x_lab, nrow = 2, rel_heights = c(1, 0.1))


# Add screenshots
sshots_plot = cowplot::ggdraw() +
  cowplot::draw_image(SSHOTS[1], x = 0, y = 1) +
  cowplot::draw_image(SSHOTS[2], x = 0.5, y = 1)

ggsave(OUT,
       final,
       device = "png",
       width = 15,
       height = 7.142,
       units = "in",
       dpi = 400)
