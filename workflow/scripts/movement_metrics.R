# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

# Debug
#IN_FILE = "/nfs/research/birney/users/ian/pilot/final_tracks/novel_object/20190615_1549_icab_hni_L_q4.csv"
#INTERVAL = 0.2
#SOURCE_FILE = here::here("workflow/scripts/movement_metrics_source.R")

## True
IN_FILE = snakemake@input[[1]]
OUT_FILE = snakemake@output[[1]]
INTERVAL = snakemake@params[["seconds_interval"]] %>% 
  as.numeric()
SOURCE_FILE = snakemake@params[["source_file"]]

source(SOURCE_FILE)

# Read in file

df = readr::read_csv(IN_FILE)

# Create vector of intervals (in seconds)

int_vec = seq(0, max(df$seconds), INTERVAL)

# Process

df = df %>% 
  tidyr::pivot_longer(cols = c(ref_x, ref_y, test_x, test_y),
                      names_to = c("fish", "axis"),
                      names_sep = "_",
                      values_to = "coord") %>% 
  # give x and y their own columns
  tidyr::pivot_wider(names_from = axis, values_from = coord) %>% 
  # split by fish
  split(., f = .$fish) %>%
  purrr::map(., function(FISH){
    out = FISH %>% 
      # ensure data frame is ordered
      dplyr::arrange(frame)
    # Get indices for rows with closest values to each value in the interval vector
    indices = purrr::map_int(int_vec, function(x){
      which.min(abs(out$seconds - x))
    })
    out = out %>% 
      # Filter df for those indices
      dplyr::slice(indices) %>% 
      # Get lagged values
      dplyr::mutate(x_lag1 = dplyr::lag(x, n = 1),
                    y_lag1 = dplyr::lag(y, n = 1),
                    x_lag2 = dplyr::lag(x, n = 2),
                    y_lag2 = dplyr::lag(y, n = 2)) %>% 
      # Get distance
      dplyr::mutate(distance = get_dist(x, x_lag1, y, y_lag1),
                    distance_b = get_dist(x, x_lag2, y, y_lag2),
                    angle = get_angle(x, x_lag1, x_lag2, y, y_lag1, y_lag2))
    
    return(out)
    }) %>% 
  # bind back into single DF
  dplyr::bind_rows()
  
# Write to file

readr::write_csv(df, OUT_FILE)


