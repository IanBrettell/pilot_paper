###################
# Libraries
###################

import os
import pandas as pd
import itertools

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

#SAMPLES_PAIRS = pd.read_csv(config["samples_pairs_file"])

## Create list of variable lists
full_list = [SAMPLES, ASSAYS, QUADRANTS]
## Create list of tuple combinations
combos = list(itertools.product(*full_list))

## Remove unavailable combinations
#excl_df = excl_df.reset_index()
#for index, row in excl_df.iterrows():
#    combos.remove((row['sample'], row['assay'], row['quadrant']))

# Create new lists of variables
SAMPLES = [i[0] for i in combos]
ASSAYS = [i[1] for i in combos]
QUADRANTS = [i[2] for i in combos]