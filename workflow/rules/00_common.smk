###################
# Libraries
###################

import os.path
import pandas as pd

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

SAMPLES_PAIRS = pd.read_csv(config["samples_pairs_file"])
