# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(magick)

# Get variables

## Debug
SETUP = here::here("book/figs/misc/setup_picture.pdf")
FRAME = here::here("book/figs/labs_and_paths_frame_grabs/0.08/dist_angle/15/open_field/20190613_1054_icab_hdr_R_110.png")
SETUP_PNG = here::here("book/figs/misc/setup_picture.png")
OUT = here::here("book/figs/misc/setup_fig/0.08/dist_angle/15/open_field/20190613_1054_icab_hdr_R/setup_pic_with_frame_110.png")


## True
SETUP = snakemake@input[["setup_pic"]]
FRAME = snakemake@input[["frame_grab"]]
SETUP_PNG = snakemake@params[["setup_png"]]
OUT = snakemake@output[["fig"]]

# Read in figures and compose

setup = magick::image_read_pdf(SETUP)
frame = magick::image_read(FRAME)

# Save setup figure as .png

magick::image_write(setup, SETUP_PNG)

# Compose image

out = cowplot::ggdraw() +
  cowplot::draw_image(SETUP_PNG,
                      x = 0, y = 0.4,
                      width = 1, height = 0.6) +
  cowplot::draw_image(FRAME,
                      x = 0, y = 0,
                      width = 1, height = 0.4) +
  # add labels
  cowplot::draw_plot_label(c("A", "B", "C"),
                           x = c(0,0,0.5),
                           y = c(1, 0.4, 0.4)) + 
  # add quadrant labels
  cowplot::draw_plot_label(c("I", "II", "III", "IV"),
                           x = c(0.46,0.22,0.22,0.46),
                           y = c(0.22,0.22,0.03,0.03),
                           color = "blueviolet",size = 8)   

## resize `frame` to make it the same width as `setup`
#
#setup_wid = 1000
#
#frame_rat = 1412/704
#
#new_hei = round(1000/frame_rat)
#
#frame_resize = magick::image_resize(frame, geometry = paste(setup_wid, new_hei, sep = "x"))
#
## Compose
#
#out = magick::image_append(c(setup, frame_resize),
#                           stack = T)
#
## Add labels
#
#out = magick::image_annotate(out, "A", location = "+10+10", style = "bold", size = 50, color = "black")
#out = magick::image_annotate(out, "B", location = "+10+760", style = "bold", size = 50, color = "black")
#out = magick::image_annotate(out, "C", location = "+510+760", style = "bold", size = 50, color = "black")
#
## Write
#
#magick::image_write(out, OUT)

ggsave(OUT,
       device = "png",
       width = 9,
       height = 11.3,
       units = "in",
       dpi = 400)
