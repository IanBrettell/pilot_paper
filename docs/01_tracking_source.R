#! /usr/bin/RScript

########################
# Libraries
########################

library(tidyverse)
library(plotly)
library(gganimate)

########################
# Functions
########################

accumulate_by <- function(dat, var) {
  var <- lazyeval::f_eval(var, dat)
  lvls <- plotly:::getLevels(var)
  dats <- lapply(seq_along(lvls), function(x) {
    cbind(dat[var %in% lvls[seq(1, x)], ], time = lvls[[x]])
  })
  dplyr::bind_rows(dats)
}

# Take list with META and TRACKING elements, and return named vector to recode fish `A` or `B` as another line
label_fishes = function(input){
  # Pull first frame from tracking data
  frame_1 = input[["TRACKING"]][1, ]
  # Pull iCab location from metadata
  ICAB_LOC = input[["META"]]$ICAB_LOC
  # create recode vector
  recode_vec = character(length = 2)
  # Set iCab as either fish A or B
  if (is.na(ICAB_LOC) == T){
    recode_vec[1] = "iCab_1"
    recode_vec[2] = "iCab_2"
  } else if (ICAB_LOC == "Top"){
    recode_vec[1] = ifelse(frame_1$Y1 < frame_1$Y2, "iCab", input[["META"]]$TEST)
    recode_vec[2] = ifelse(frame_1$Y2 < frame_1$Y1, "iCab", input[["META"]]$TEST)
  } else if (ICAB_LOC == "Bottom"){
    recode_vec[1] = ifelse(frame_1$Y1 > frame_1$Y2, "iCab", input[["META"]]$TEST)
    recode_vec[2] = ifelse(frame_1$Y2 > frame_1$Y1, "iCab", input[["META"]]$TEST)
  } else if (ICAB_LOC == "Left"){
    recode_vec[1] = ifelse(frame_1$X1 < frame_1$X2, "iCab", input[["META"]]$TEST)
    recode_vec[2] = ifelse(frame_1$X2 < frame_1$X1, "iCab", input[["META"]]$TEST)
  } else if (ICAB_LOC == "Right"){
    recode_vec[1] = ifelse(frame_1$X1 > frame_1$X2, "iCab", input[["META"]]$TEST)
    recode_vec[2] = ifelse(frame_1$X2 > frame_1$X1, "iCab", input[["META"]]$TEST)
  }
  
  # Recode the recode vector to get right capitalisation of line names
  line_names = c("iCab", "HdrR", "Kaga", "HNI", "HO5", "iCab_1", "iCab_2")
  names(line_names) = c("iCab", "hdr", "kaga", "hni", "ho5", "iCab_1", "iCab_2")
  recode_vec = dplyr::recode(recode_vec, !!!line_names)
  # Set names and return vector
  names(recode_vec) = c("A", "B")
  return(recode_vec)
}

# Generate tracking plot

tracking_plot = function(input){
  # Create long DF 
  df = input %>%
    accumulate_by(~TIME)
  # Generate plot
  out_plot = df %>% 
    plot_ly(x = ~COORD_X, y = ~COORD_Y) %>% 
    add_paths(color = ~LINE,
              colors = fish_pal,
              line = list(simplify = F),
              frame = ~time) %>% 
    layout(yaxis = list(autorange = "reversed")) %>% 
    layout(xaxis = list(scaleanchor = "y", scaleratio = 1))
  
  return(out_plot)
}

########################
# Plotting parameters
########################

# `darker` and `lighter` from `karyoploteR`

source("https://gist.githubusercontent.com/brettellebi/c5015ee666cdf8d9f7e25fa3c8063c99/raw/91e601f82da6c614b4983d8afc4ef399fa58ed4b/karyoploteR_lighter_darker.R")

# Palettes

## iCab control colours ascertained with:
## scales::show_col(c(darker(fish_pal[1], amount = 50), lighter(fish_pal[1], amount = 25)))

fish_pal = c("#F1BB7B", "#FA796C", "#AC3E3F", "#79301F", "#D67236", "#BF8949", "#FFD494")
names(fish_pal) = c("iCab", "HdrR", "Kaga", "HNI", "HO5", "iCab_1", "iCab_2")
