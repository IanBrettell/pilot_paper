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
import itertools
from collections import OrderedDict

# Get variables

VID = "/hps/nobackup/birney/users/ian/pilot/stiched/open_field/20190613_1307_icab_hdr_R.avi"
PATH = "/hps/nobackup/birney/users/ian/pilot/tmp_path.avi"
HMM = [
    "/hps/nobackup/birney/users/ian/pilot/hmm_path_videos/0.08/dist_angle/15/open_field/20190613_1307_icab_hdr_R_ref.avi",
    "/hps/nobackup/birney/users/ian/pilot/hmm_path_videos/0.08/dist_angle/15/open_field/20190613_1307_icab_hdr_R_test.avi"
]
HMM.sort()

# Capture videos

cap_vid = cv2.VideoCapture(VID)
cap_path = cv2.VideoCapture(PATH)
cap_hmm_ref = cv2.VideoCapture(HMM[0])
cap_hmm_test = cv2.VideoCapture(HMM[1])

####################
# Get N frames of video, width, and height
####################

test = {"vid" : {'n_frames' : 18000, 'wid' : 704, 'hei' : 706}}

IN = [VID, PATH, HMM[0], HMM[1]]
IN = {
    "vid": {"path" : VID},
    "path": {"path": PATH},
    "hmm_ref": {"path" : HMM[0]},
    "hmm_test": {"path" : HMM[1]}
    }

for key in IN:
    cap = cv2.VideoCapture(IN[key]['path'])
    IN[key]['n_frames'] = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    IN[key]['fps'] = int(cap.get(cv2.CAP_PROP_FPS))
    IN[key]['wid'] = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    IN[key]['hei'] = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))

####################
# Downsample frames for path, hmm_ref, and hmm_test videos
####################

ind_dict = {
    'path': {},
    'hmm_ref': {},
    'hmm_test': {}
}

for key in ind_dict:
    ratio = round(IN[key]['n_frames']/IN['vid']['n_frames'])
    keep_old = list(range(0,IN[key]['n_frames']))[::ratio]
    keep_new = keep_old
    extra = IN['vid']['n_frames'] - len(keep_old)
    elements = list(itertools.repeat(keep_old[len(keep_old)-1], extra))
    for i in elements:
        keep_new.append(i)
    ind_dict[key]['ind_to_keep'] = keep_new


####################
# Resize
####################

cap_path.set(1, 0)
ret_path, frame_path = cap_path.read()
test = cv2.resize(frame_path, (WID_VID, HEI_VID))

####################
# Set up video writer
####################

fourcc = cv2.VideoWriter_fourcc('h', '2', '6', '4')

## Debug
video_writer = cv2.VideoWriter(
    "/hps/nobackup/birney/users/ian/pilot/tmp_out.avi",
    fourcc,
    IN['vid']['fps'],
    (IN['vid']['wid']*2, IN['vid']['hei']),
    isColor = True
)

## True
video_writer = cv2.VideoWriter(
    OUT,
    fourcc,
    IN['vid']['fps'],
    (IN['vid']['wid']*2, IN['vid']['hei']*2),
    isColor = True
)

####################
# Resize and write
####################

start = 0
end = 30
i = start
while i in range(start,end):
    cap_vid.set(1, i)
    ret_vid, frame_vid = cap_vid.read()
    target_frame = ind_dict['path']['ind_to_keep'][i]
    cap_path.set(1, target_frame)
    ret_path, frame_path = cap_path.read()
    frame_path_new = cv2.resize(frame_path, (IN['vid']['wid'], IN['vid']['hei']))
    top = np.hstack((frame_vid, frame_path_new))
    i += 1

while i in range(start,end):
    cap_vid.set(1, i)
    ret_vid, frame_vid = cap_vid.read()
    target_frame = ind_dict['path']['ind_to_keep'][i]
    cap_path.set(1, target_frame)
    ret_path, frame_path = cap_path.read()
    frame_path_new = cv2.resize(frame_path, (IN['vid']['wid'], IN['vid']['hei']))
    top = np.hstack((frame_vid, frame_path_new))
    video_writer.write(top)
    i += 1