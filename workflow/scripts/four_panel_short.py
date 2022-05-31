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

# Import libraries

import numpy as np
import pandas as pd
import cv2
import os
import imageio
import pygifsicle

# Get variables

## Debug
IN = "/hps/nobackup/birney/users/ian/pilot/four_panel_vids/0.08/dist_angle/15/open_field/20190613_1338_icab_ho5_R.avi"
TOT_SEC = 60
OUT = "/hps/nobackup/birney/users/ian/pilot/tmp_out.avi"

# True
IN = snakemake.input[0]
TOT_SEC = int(snakemake.params.tot_sec)
OUT = snakemake.output.avi

########################
## Get video params
########################

# Capture video

cap = cv2.VideoCapture(IN)

# Number of frames

N_FRAMES = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
FPS = int(cap.get(cv2.CAP_PROP_FPS))
WID = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
HEI = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))

# Get frames to write

FPS_NEW = TOT_SEC * FPS
SPLIT = int(N_FRAMES/FPS_NEW)
T_FRAMES = list(range(N_FRAMES))[::SPLIT]


#######################
# Video writer
#######################

fourcc = cv2.VideoWriter_fourcc('h', '2', '6', '4')

video_writer = cv2.VideoWriter(
    OUT,
    fourcc,
    FPS,
    (WID, HEI),
    isColor = True
)

#######################
# Write video
#######################

for i in T_FRAMES:
    # Set frame
    cap.set(cv2.CAP_PROP_POS_FRAMES, i)
    # Read frame
    ret_vid, frame = cap.read()
    # Write to video
    video_writer.write(frame)

video_writer.release()

########################
## Add frames to gif
########################
#
#out_list = []
#for i in T_FRAMES:
#    # Set frame
#    cap.set(cv2.CAP_PROP_POS_FRAMES, i)
#    # Read frame
#    ret_vid, frame = cap.read()
#    # Add to list
#    out_list.append(i)
#
########################
## Write to gif
########################
#
#imageio.mimsave(OUT, out_list, fps = 3)
#pygifsicle.optimize(OUT)

