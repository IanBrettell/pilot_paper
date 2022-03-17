###################
# Libraries
###################

import os
import pandas as pd
import itertools
import numpy as np

###################
# Config
###################

configfile: "config/config.yaml"

###################
# Variables
###################

samples_df = pd.read_csv(config["samples_file"], comment = '#')

SAMPLES = samples_df["sample"]
ASSAYS = ["open_field", "novel_object"]
QUADRANTS = ["q1", "q2", "q3", "q4"]

# Read in videos to be excluded
excl_df = pd.read_csv(config["excluded_videos"], comment = "#")

## Create list of variable lists
full_list = [SAMPLES, ASSAYS, QUADRANTS]
## Create list of tuple combinations
combos = list(itertools.product(*full_list))

# Remove unavailable combinations
excl_df = excl_df.reset_index()
for index, row in excl_df.iterrows():
    combos.remove((row['sample'], row['assay'], row['quadrant']))

# Create new lists of variables
SAMPLES = [i[0] for i in combos]
ASSAYS = [i[1] for i in combos]
QUADRANTS = [i[2] for i in combos]

# Multiply lists for combinations with `seconds_interval`
n_intervals = len(config["seconds_interval"])
SAMPLES_INT = SAMPLES * n_intervals
ASSAYS_INT = ASSAYS * n_intervals
QUADRANTS_INT = QUADRANTS * n_intervals
## Multiply each element of `seconds_intervals` by length of original SAMPLES list
INTERVALS_INT = np.repeat(config["seconds_interval"], len(SAMPLES))