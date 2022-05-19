# Send log
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
import cv2

# Get variables

## Debug

#IN = "/hps/nobackup/birney/users/ian/pilot/path_frames/0.08/dist_angle/15/open_field/20190612_1236_icab_kaga_R/18000.png"
#DIMS = "/hps/software/users/birney/ian/repos/pilot_paper/config/split_video_dims.csv"
#FPS = 30
#ASSAY = "open_field"
#SAMPLE = "20190612_1236_icab_kaga_R"

## True

IN = snakemake.input.frame_path
DIMS = snakemake.input.dims[0]
FPS = int(snakemake.params.fps)
ASSAY = snakemake.params.assay
SAMPLE = snakemake.params.sample
OUT = snakemake.output[0]

# Read in dims file

dims = pd.read_csv(DIMS)

## Filter for target assay and sample (creates 4 rows, 1 per quadrant)

dims_filt = dims.loc[(dims["assay"] == ASSAY) & (dims["sample"] == SAMPLE)]

## Get N frames

N_FRAMES = max(dims_filt['n_frames'])

## Get video width and height

TOT_WID = dims_filt.loc[(dims_filt["quadrant"] == "q1") | (dims_filt["quadrant"] == "q2")]['wid'].sum()
TOT_HEI = dims_filt.loc[(dims_filt["quadrant"] == "q1") | (dims_filt["quadrant"] == "q4")]['hei'].sum()

# Get list of files

DIRNAME = os.path.dirname(IN)
frame_files = []
for i in range(1,N_FRAMES + 1):
    frame_files.append(
        os.path.join(DIRNAME, str(i) + ".png")
    )

####################
# Set up video writer
####################

fourcc = cv2.VideoWriter_fourcc('h', '2', '6', '4')

video_writer = cv2.VideoWriter(
    OUT,
    fourcc,
    FPS,
    (TOT_WID, TOT_HEI),
    isColor = True
)

####################
# Write video
####################

start = 0
end = N_FRAMES 
i = start
while i in range(start, end):
    frame = cv2.imread(frame_files[i])
    # resize
    out = cv2.resize(frame, (TOT_WID, TOT_HEI))
    video_writer.write(out)
    i += 1

video_writer.release()
