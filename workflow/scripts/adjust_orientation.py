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
SAMPLES_FILE = "config/samples.csv"
OUT_TEST = "results/test.png"
OUT_FILE = "/hps/nobackup/birney/users/ian/pilot/flipped/20190611_1410_icab_kaga_R.avi"

## True
IN_FILE = snakemake.input[0]
SAMPLE = snakemake.params.sample
SAMPLES_FILE = snakemake.params.samples_file
OUT_FILE = snakemake.output[0]

# Read samples_file

samples_df = pd.read_csv(SAMPLES_FILE, comment="#", skip_blank_lines=True, index_col=0)

# Get tank side

tank_side = samples_df.loc[SAMPLE, "tank_side"]

# Read video from file

cap = cv.VideoCapture(IN_FILE)

# Frame width and height

wid = int(cap.get(cv.CAP_PROP_FRAME_WIDTH))
hei = int(cap.get(cv.CAP_PROP_FRAME_HEIGHT))
size = (wid, hei)

# Get total frame length

vid_len = int(cap.get(cv.CAP_PROP_FRAME_COUNT))
start = 0
end = vid_len

# Get frames per second

fps = int(cap.get(cv.CAP_PROP_FPS))

# Define the codec and create VideoWriter object

fourcc = cv.VideoWriter_fourcc('h', '2', '6', '4')
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
    # Flip if tank side is "R"
    if tank_side == 'R':
        frame = cv.rotate(frame, cv.ROTATE_180)
    # Write frame
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


