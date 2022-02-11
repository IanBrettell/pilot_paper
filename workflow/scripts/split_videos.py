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

# Read samples_file

samples_df = pd.read_csv(snakemake.params.samples_file, comment="#", skip_blank_lines=True, index_col=0)

## Get start and end frames

if snakemake.wildcards.assay == "open_field":
    start = int(samples_df.loc[snakemake.wildcards.sample, "of_start"])
    end = int(samples_df.loc[snakemake.wildcards.sample, "of_end"])
elif snakemake.wildcards.assay == "novel_object":
    start = int(samples_df.loc[snakemake.wildcards.sample, "no_start"])
    end = int(samples_df.loc[snakemake.wildcards.sample, "no_end"])

# Get crop adjustment values

adj_top = int(samples_df.loc[snakemake.wildcards.sample, "adj_top"])
adj_bottom = int(samples_df.loc[snakemake.wildcards.sample, "adj_bottom"])
adj_left = int(samples_df.loc[snakemake.wildcards.sample, "adj_left"])
adj_right = int(samples_df.loc[snakemake.wildcards.sample, "adj_right"])

# Get key variables
in_file = snakemake.input[0]
#start = int(snakemake.params.start)
#end = int(snakemake.params.end)
quadrant = snakemake.wildcards.quadrant
out_file = snakemake.output[0]

# Read video from file
cap = cv.VideoCapture(in_file)

# Frame width and height
wid = int(cap.get(cv.CAP_PROP_FRAME_WIDTH))
hei = int(cap.get(cv.CAP_PROP_FRAME_HEIGHT))
# Get total frame length
vid_len = int(cap.get(cv.CAP_PROP_FRAME_COUNT))
# Get frames per second
fps = int(cap.get(cv.CAP_PROP_FPS))
# Adapt basename prefix
in_file_pref = os.path.split(in_file)[1].split('.')[0]
# Get file name without extension
meta_list = os.path.split(in_file)[1].split('.')[0].split('_')

# Get metadata
date = meta_list[0]
time = meta_list[1]
test_line = meta_list[3]
tank_side = meta_list[4]

# Get bounding box coords for target quadrant
if quadrant == 'q1':
    top = 0
    bottom = round(((hei - 1) / 2) + adj_bottom)
    left = round(((wid - 1) / 2) + adj_left)
    right = wid - 1
elif quadrant == 'q2':
    top = 0
    bottom = round(((hei - 1) / 2) + adj_bottom)
    left = 0
    right = round(((wid - 1) / 2) + adj_right)
elif quadrant == 'q3':
    top = round(((hei - 1) / 2) + adj_top)
    bottom = hei - 1
    left = 0
    right = round(((wid - 1) / 2) + adj_right)
elif quadrant == 'q4':
    top = round(((hei - 1) / 2) + adj_top)
    bottom = hei - 1
    left = round(((wid - 1) / 2) + adj_left)
    right = wid  - 1
else:
    print('Invalid quadrant')

# Get size of output video
size = (right - left, bottom - top)

# Define the codec and create VideoWriter object
# Works on local Mac:
fourcc = cv.VideoWriter_fourcc('m', 'p', '4', 'v')
out = cv.VideoWriter(out_file, fourcc, fps, size, isColor=True)

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
    # Flip if tank side is "R"
    if tank_side == 'R':
        frame = cv.rotate(frame, cv.ROTATE_180)
    # Crop frame
    frame = frame[top:bottom, left:right]
    # Write frame
    out.write(frame)
    # Show image
#    cv.imshow('frame', frame)
    # Add to counter
    i += 1
    # Press 'esc' to close video
#    if cv.waitKey(1) == 27:
#        cv.destroyAllWindows()
#        cv.waitKey(1)
#        break

cap.release()
out.release()
out = None
#cv.destroyAllWindows()
#cv.waitKey(1)
