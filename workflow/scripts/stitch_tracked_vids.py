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

import cv2 as cv
import numpy as np

# Get variables

## Debug

IN = ["/hps/nobackup/birney/users/ian/pilot/tracked/novel_object/20190616_1045_icab_icab_L_q1.avi",
      "/hps/nobackup/birney/users/ian/pilot/tracked/novel_object/20190616_1045_icab_icab_L_q3.avi",
      "/hps/nobackup/birney/users/ian/pilot/tracked/novel_object/20190616_1045_icab_icab_L_q2.avi",
      "/hps/nobackup/birney/users/ian/pilot/tracked/novel_object/20190616_1045_icab_icab_L_q4.avi"
]
#FPS = 30
#OUT = "/hps/nobackup/birney/users/ian/pilot/tmp.avi"

## True
IN = snakemake.input
OUT = snakemake.output[0]
FPS = int(snakemake.params.fps)

# Sort order to ensure q1 is first, etc.

IN = sorted(IN)

# Capture videos

cap1 = cv.VideoCapture(IN[0])
cap2 = cv.VideoCapture(IN[1])
cap3 = cv.VideoCapture(IN[2])
cap4 = cv.VideoCapture(IN[3])

## Set up fourcc

fourcc = cv.VideoWriter_fourcc('h', '2', '6', '4')

# Get width and height

WID = int(cap1.get(cv.CAP_PROP_FRAME_WIDTH) + cap2.get(cv.CAP_PROP_FRAME_WIDTH))
HEI = int(cap1.get(cv.CAP_PROP_FRAME_HEIGHT) + cap4.get(cv.CAP_PROP_FRAME_HEIGHT))

# Get total video length

VID_LEN = int(cap1.get(cv.CAP_PROP_FRAME_COUNT))

# Set up video writer
video_writer = cv.VideoWriter(
    OUT,
    fourcc,
    FPS,
    (WID, HEI),
    isColor = True
)

# Compile all four videos into one 

start = 0
end = VID_LEN
i = start
while i in range(start,end):
    cap1.set(1, i)
    ret1, frame1 = cap1.read()
    cap2.set(1, i)
    ret2, frame2 = cap2.read()
    cap3.set(1, i)
    ret3, frame3 = cap3.read()
    cap4.set(1, i)
    ret4, frame4 = cap4.read()
    top = np.hstack((frame2, frame1))
    bottom = np.hstack((frame3, frame4))
    out = np.vstack((top, bottom))
    video_writer.write(out)
    i += 1

# Release captures

cap1.release()
cap2.release()
cap3.release()
cap4.release()

# Release video

video_writer.release()
