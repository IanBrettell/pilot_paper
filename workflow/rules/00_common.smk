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

SAMPLES_ALL = samples_df["sample"]
ASSAYS_ALL = ["open_field", "novel_object"]
QUADRANTS_ALL = ["q1", "q2", "q3", "q4"]
REF_TEST = ["ref", "test"]

# Read in videos to be excluded
excl_df = pd.read_csv(config["excluded_videos"], comment = "#")

## Create list of variable lists
full_list = [SAMPLES_ALL, ASSAYS_ALL, QUADRANTS_ALL]
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

# Create combinations of assay and sample, excluding 20190616_1227_icab_kaga_R
assay_sample_list = [SAMPLES_ALL, ASSAYS_ALL]
assay_sample_combos = list(itertools.product(*assay_sample_list))
for index, row in excl_df.iterrows():
    assay_sample_combos.remove((row['sample'], row['assay']))
AS_SAMPLES = [i[0] for i in assay_sample_combos]
AS_ASSAYS = [i[1] for i in assay_sample_combos]