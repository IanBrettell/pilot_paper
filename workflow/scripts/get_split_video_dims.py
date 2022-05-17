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

# Get variables

## Debug
IN = ["/hps/nobackup/birney/users/ian/pilot/split/open_field/20190611_1331_icab_icab_R_q1.avi",
      "/hps/nobackup/birney/users/ian/pilot/split/open_field/20190611_1331_icab_icab_R_q2.avi",
      "/hps/nobackup/birney/users/ian/pilot/split/open_field/20190611_1331_icab_icab_R_q3.avi"]

## True
IN = snakemake.input
OUT = snakemake.output[0]

# Sort file names

IN = sorted(IN)

# Write lines to file and close

file = open(OUT, 'w')

# write header
file.writelines("assay,sample,quadrant,wid,hei,n_frames\n")

for VID in IN:
    cap = cv.VideoCapture(VID)
    wid = str(int(cap.get(cv.CAP_PROP_FRAME_WIDTH)))
    hei = str(int(cap.get(cv.CAP_PROP_FRAME_HEIGHT)))
    n_frames = str(int(cap.get(cv.CAP_PROP_FRAME_COUNT)))
    # get sample
    sample = '_'.join(os.path.basename(VID).replace('.avi', '').split('_')[:-1])
    # get quadrant
    quadrant = os.path.basename(VID).replace('.avi', '').split('_')[-1]
    # get assay
    assay = os.path.dirname(VID).split('/')[-1]
    # compose full line
    line_to_write = ','.join([assay,sample,quadrant,wid,hei,n_frames]) + '\n'
    file.writelines(line_to_write)

file.close()
