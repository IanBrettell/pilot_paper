#!/usr/bin/env python3

# example bash script: python ~/Documents/Repositories/pilot_paper/code/scripts/20201203_amend_h5s.py /Users/brettell/Desktop/labeled-data/20190616_1717_icab_icab_L_q4_novel_object/CollectedData_ian_brettell.h5

# Import libraries
import sys
import os
import pandas as pd

# Get file name
in_file = sys.argv[1]

# Read in h5 file
df = pd.read_hdf(in_file, "df_with_missing")

# Find columns with all missing
test_lines = ['icab_b', 'hni', 'kaga', 'hdrr', 'ho5', 'mikk_a', 'mikk_b']

lines_to_drop = []
for line in test_lines:
    if pd.isnull(df['ian_brettell'][line]).all().all():
        lines_to_drop.append(line)

# Drop them from DF
df.drop(columns = lines_to_drop, level = 1, inplace = True)

# Rename columns
pairs = [['icab_a', 'ref'], ['icab_b', 'test'], ['hni', 'test'], ['kaga', 'test'], ['hdrr', 'test'], ['ho5', 'test'], ['mikk_a', 'test'], ['mikk_b', 'test']]
dict_ = {k: v for k, v in pairs}
df.rename(dict_, level="individuals", axis=1, inplace=True)

# Write file
os.remove(in_file) # remove old file
df.to_hdf(in_file, key = "df_with_missing")
