###################
# Libraries
###################

import os.path
import pandas as pd

###################
# Config
###################

configfile: "code/snakemake/20210701/config/config.yaml"

###################
# Variables
###################

SAMPLES = pd.read_csv(config["samples_file"], comment="#", skip_blank_lines=True, index_col=0)
ASSAYS = ["open_field", "novel_object"]
QUADRANTS = ["q1", "q2", "q3", "q4"]

