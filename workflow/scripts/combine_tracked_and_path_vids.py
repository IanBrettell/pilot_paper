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

import numpy as np
import cv2

# Get variables

HMM = [
    "/hps/nobackup/birney/users/ian/pilot/hmm_path_videos/0.08/dist_angle/15/open_field/20190613_1307_icab_hdr_R_ref.avi",
    "/hps/nobackup/birney/users/ian/pilot/hmm_path_videos/0.08/dist_angle/15/open_field/20190613_1307_icab_hdr_R_test.avi"
]
