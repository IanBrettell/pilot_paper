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

# Get variables

## Debug
LABELS = "/hps/nobackup/birney/users/ian/pilot/stitched/open_field/20190613_1054_icab_hdr_R.avi"
PATHS = "/hps/nobackup/birney/users/ian/pilot/path_vids/0.08/dist_angle/15/open_field/20190613_1054_icab_hdr_R.avi"
SECOND = 110
OUT = "/hps/nobackup/birney/users/ian/pilot/tmp_out.png"

# True
LABELS = snakemake.input.labels
PATHS = snakemake.input.paths
SECOND = int(snakemake.params.target_second)
OUT = snakemake.output[0]

########################
## Capture videos
########################

cap_vid = cv2.VideoCapture(LABELS)
cap_path = cv2.VideoCapture(PATHS)

####################
# Get N frames of video, width, and height
####################

IN = {
    "vid": {"path" : LABELS},
    "path": {"path": PATHS},
    }

for key in IN:
    cap = cv2.VideoCapture(IN[key]['path'])
    IN[key]['n_frames'] = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    IN[key]['fps'] = int(cap.get(cv2.CAP_PROP_FPS))
    IN[key]['wid'] = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    IN[key]['hei'] = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))


# Throw error if n_frames aren't equal for all videos

if not (IN['vid']['n_frames'] == IN['path']['n_frames']):
    raise Exception('Videos do not have the same number of frames')

# Get total height and width in pixels (by doubling original video)

TOT_WID = IN['vid']['wid']*2
TOT_HEI = IN['vid']['hei']

# Get FPS

FPS = IN['vid']['fps']

# Get target frame

FRAME = FPS*SECOND

#######################
# Write frame
#######################

# Set frames
cap_vid.set(cv2.CAP_PROP_POS_FRAMES, FRAME)
cap_path.set(cv2.CAP_PROP_POS_FRAMES, FRAME)

# Capture frames
ret_vid, frame_vid = cap_vid.read()
ret_path, frame_path = cap_path.read()

# Stack frames horizontally
out = np.hstack((frame_vid, frame_path))

# Write to file
cv2.imwrite(OUT, out)
