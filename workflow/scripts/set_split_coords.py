# Send stdout and stderr to log file
import sys,os
import logging, traceback
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    )
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                         *traceback.format_exception(exc_type, exc_value, exc_traceback)
                         ])
                 )
# Install exception handler
sys.excepthook = handle_exception

# Import libraries
import pandas as pd
import numpy as np
import cv2 as cv
import os

# Get variables

## Debugging
IN_FILE = "/nfs/research/birney/projects/indigene/raw_data/ian_videos/ian_pilot/all/20190612_1053_icab_kaga_R.avi"
SAMPLES_FILE = "config/samples.csv"
SAMPLE = "20190612_1053_icab_kaga_R"
OUT_FILE = "results/split_coord_images/20190612_1053_icab_kaga_R.png"

## True
IN_FILE = snakemake.input.video
SAMPLES_FILE = snakemake.params.samples_file
SAMPLE = snakemake.params.sample
OUT_FILE = snakemake.output.fig

# Read samples_file
samples_df = pd.read_csv(SAMPLES_FILE, comment="#", skip_blank_lines=True, index_col=0)

# Get date
date = int(samples_df.loc[SAMPLE, "date"])

# Get start frame for open field assay
start = int(samples_df.loc[SAMPLE, "of_start"])

# Get crop adjustment values
## note: Negative values for top/bottom shift boundary up
## note: Negative values for left/right shift boundary left
adj_top = int(samples_df.loc[SAMPLE, "adj_top"])
adj_bottom = int(samples_df.loc[SAMPLE, "adj_bottom"])
adj_left = int(samples_df.loc[SAMPLE, "adj_left"])
adj_right = int(samples_df.loc[SAMPLE, "adj_right"])

# Read video from file
cap = cv.VideoCapture(IN_FILE)

# Frame width and height
wid = int(cap.get(cv.CAP_PROP_FRAME_WIDTH))
hei = int(cap.get(cv.CAP_PROP_FRAME_HEIGHT))

# Set adjusted midpoints
mid_x = round(((wid - 1) / 2) + adj_right)
mid_y = round(((hei - 1) / 2) + adj_top)

# Capture start frame
cap.set(cv.CAP_PROP_POS_FRAMES, start)

# Read frame
ret, frame = cap.read()

# Add vertical line 
start_point = (mid_x, 0)
end_point = (mid_x, hei)
color = (255,0,0)
thickness = 1
frame = cv.line(frame, start_point, end_point, color, thickness)

# Add horizontal line
start_point = (0, mid_y)
end_point = (wid, mid_y)
color = (255,0,0)
thickness = 1
frame = cv.line(frame, start_point, end_point, color, thickness)

# Write frame
cv.imwrite(OUT_FILE, frame)

cap.release()
