# Model covariates


```r
library(tidyverse)
library(ggbeeswarm)
library(devtools)
library(cowplot)
library(ggpubr)
```

## Input file paths


```r
IN = "/hps/nobackup/birney/users/ian/pilot/hmm_out/0.05/dist_angle/18.csv"

# Get lighter/darker functions

devtools::source_gist("c5015ee666cdf8d9f7e25fa3c8063c99")
#> ℹ Sourcing https://gist.githubusercontent.com/brettellebi/c5015ee666cdf8d9f7e25fa3c8063c99/raw/15832e2684e4c08a652eb82d4b559bea4e8994e4/karyoploteR_lighter_darker.R
#> ℹ SHA-1 hash of file is 3d840ff655bd54f4b19517cf99ec4bd3640a30c5
```

## Read in file


```r
# Read in file

df = readr::read_csv(IN) 
#> Rows: 14649781 Columns: 15
#> ── Column specification ────────────────────────────────────
#> Delimiter: ","
#> chr (6): assay, ref_fish, test_fish, tank_side, quadrant...
#> dbl (9): date, time, frame, seconds, x, y, distance, ang...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

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
#> `summarise()` has grouped output by 'assay', 'date',
#> 'time', 'quadrant', 'tank_side'. You can override using the
#> `.groups` argument.

# Get mean speed over both assays for each individual

df_control_noassay = df %>% 
  # Filter for only iCabs paired with iCabs
  dplyr::filter(test_fish == "icab") %>% 
  # Get individual
  tidyr::unite(date, time, quadrant, fish,
               col = "indiv", 
               remove = F) %>% 
  # Group by assay and individual to get mean speed
  dplyr::group_by(date, time, quadrant, tank_side, indiv) %>% 
  # Calculate mean speed
  dplyr::summarise(mean_speed = mean(distance)) %>% 
  dplyr::ungroup() %>% 
  # Make date a factor
  dplyr::mutate(date = as.factor(date))
#> `summarise()` has grouped output by 'date', 'time',
#> 'quadrant', 'tank_side'. You can override using the
#> `.groups` argument.

###
```


