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
import sys

# Get variables

## Debugging
IN_FILE = "/nfs/research/birney/projects/indigene/raw_data/ian_videos/ian_pilot/all/20190611_1410_icab_kaga_R.avi"
SAMPLE = "20190611_1410_icab_kaga_R"
ASSAY = "novel_object"
QUADRANT = "q1"
SAMPLES_FILE = "config/samples.csv"
OUT_FILE = " /nfs/research/birney/users/ian/pilot/split/novel_object/20190611_1410_icab_kaga_R_q1.mp4"

## True
IN_FILE = snakemake.input[0]
SAMPLE = snakemake.params.sample
ASSAY = snakemake.params.assay
QUADRANT = snakemake.params.quadrant
SAMPLES_FILE = snakemake.params.samples_file
OUT_FILE = snakemake.output[0]

# Read samples_file

samples_df = pd.read_csv(SAMPLES_FILE, comment="#", skip_blank_lines=True, index_col=0)

# Get start and end frames

if ASSAY == "open_field":
    start = int(samples_df.loc[SAMPLE, "of_start"])
    end = int(samples_df.loc[SAMPLE, "of_end"])
elif ASSAY == "novel_object":
    start = int(samples_df.loc[SAMPLE, "no_start"])
    end = int(samples_df.loc[SAMPLE, "no_end"])

# Get crop adjustment values

adj_top = int(samples_df.loc[SAMPLE, "adj_top"])
adj_right = int(samples_df.loc[SAMPLE, "adj_right"])

# Get tank side

tank_side = samples_df.loc[SAMPLE, "tank_side"]

# Read video from file

cap = cv.VideoCapture(IN_FILE)

# Frame width and height

wid = int(cap.get(cv.CAP_PROP_FRAME_WIDTH))
hei = int(cap.get(cv.CAP_PROP_FRAME_HEIGHT))

# Get total frame length

vid_len = int(cap.get(cv.CAP_PROP_FRAME_COUNT))
# Get frames per second

fps = int(cap.get(cv.CAP_PROP_FPS))

# Set adjusted midpoints
mid_x = round(((wid - 1) / 2) + adj_right)
mid_y = round(((hei - 1) / 2) + adj_top)

# Get bounding box coords for target quadrant
if tank_side == "L":
    if QUADRANT == 'q1':
        top = 0
        bottom = mid_y
        left = mid_x
        right = wid - 1
    elif QUADRANT == 'q2':
        top = 0
        bottom = mid_y
        left = 0
        right = mid_x
    elif QUADRANT == 'q3':
        top = mid_y
        bottom = hei - 1
        left = 0
        right = mid_x
    elif QUADRANT == 'q4':
        top = mid_y
        bottom = hei - 1
        left = mid_x
        right = wid  - 1
    else:
        print('Invalid quadrant')
# NOTE:  when tank_side == "R" the video is upside-down, so q1 = q3 and q2 = q4, and vice versa.
# We then flip the video 180 degrees when extracting the video below.
elif tank_side == "R":
    if QUADRANT == 'q1':
        top = mid_y
        bottom = hei - 1
        left = 0
        right = mid_x
    elif QUADRANT == 'q2':
        top = mid_y
        bottom = hei - 1
        left = mid_x
        right = wid  - 1
    elif QUADRANT == 'q3':
        top = 0
        bottom = mid_y
        left = mid_x
        right = wid - 1
    elif QUADRANT == 'q4':
        top = 0
        bottom = mid_y
        left = 0
        right = mid_x
    else:
        print('Invalid quadrant')    
    
# Get size of output video

size = (right - left, bottom - top)

# Define the codec and create VideoWriter object

fourcc = cv.VideoWriter_fourcc('m', 'p', '4', 'v')
out = cv.VideoWriter(OUT_FILE, fourcc, fps, size, isColor=True)

# Capture frame-by-frame

i = start
while i in range(start,end):
    cap.set(cv.CAP_PROP_POS_FRAMES, i)
    # Capture frame-by-frame
    ret, frame = cap.read()
    # if frame is read correctly ret is True
    if not ret:
        print("Can't receive frame (stream end?). Exiting ...")
        break
    # Crop frame
    frame = frame[top:bottom, left:right]
    # Write frame

    # Flip if tank side is "R"
    if tank_side == 'R':
        frame = cv.rotate(frame, cv.ROTATE_180)

    out.write(frame)
    # Add to counter
    i += 1
    # Press 'esc' to close video
#    if cv.waitKey(1) == 27:
#        cv.destroyAllWindows()
#        cv.waitKey(1)
#        break

cap.release()
out.release()
#out = None
#cv.destroyAllWindows()
#cv.waitKey(1)


