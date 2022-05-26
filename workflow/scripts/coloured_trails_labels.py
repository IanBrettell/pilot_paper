###############################################################################
# This is an adaptation of the script trajectories_to_video.py part
# of the idtrackerai package.
#
# It was slightly modified in order to make it play nicely in my snakemake
# workflow
#
# Author of the modification: Saul Pierotti, saul@ebi.ac.uk
# Edited: February 2nd 2022
###############################################################################

# This file is part of idtracker.ai a multiple animals tracking system
# described in [1].
# Copyright (C) 2017- Francisco Romero Ferrero, Mattia G. Bergomi,
# Francisco J.H. Heras, Robert Hinz, Gonzalo G. de Polavieja and the
# Champalimaud Foundation.
#
# idtracker.ai is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details. In addition, we require
# derivatives or applications to acknowledge the authors by citing [1].
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# For more information please send an email (idtrackerai@gmail.com) or
# use the tools available at https://gitlab.com/polavieja_lab/idtrackerai.git.
#
# [1] Romero-Ferrero, F., Bergomi, M.G., Hinz, R.C., Heras, F.J.H., de Polavieja, G.G., Nature Methods, 2019.
# idtracker.ai: tracking all individuals in small or large collectives of unmarked animals.
# (F.R.-F. and M.G.B. contributed equally to this work.
# Correspondence should be addressed to G.G.d.P: gonzalo.polavieja@neuro.fchampalimaud.org)

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
from tqdm import tqdm
from idtrackerai.utils.py_utils import get_spaced_colors_util

## Debug
#video_object_path = "/hps/nobackup/birney/users/ian/pilot/split/open_field/session_20190611_1331_icab_icab_R_q1/video_object.npy"
#trajectories_path = "/hps/nobackup/birney/users/ian/pilot/split/open_field/session_20190611_1331_icab_icab_R_q1/trajectories_wo_gaps/trajectories_wo_gaps.npy"
#path_to_save_video = "/hps/nobackup/birney/users/ian/pilot/tmp.avi"
#SAMPLE = "20190611_1331_icab_icab_R"
#REF_LOC = "NA"

## True
video_object_path = snakemake.input.video_object
trajectories_path = snakemake.input.trajectories
path_to_save_video = snakemake.output[0]
SAMPLE = snakemake.params.sample
REF_LOC = snakemake.params.ref_loc


# Get identity of test fish
TEST_FISH = SAMPLE.split('_')[3]

# Recode line
recode_dict = {
    "icab" : "iCab",
    "hdr" : "HdrR",
    "hni" : "HNI",
    "kaga" : "Kaga",
    "ho5" : "HO5"
}

TEST_FISH = recode_dict[TEST_FISH]

# From `main()`

########################
# Read in video object and trajectories
########################

## Load video object
video_object = np.load(
        video_object_path, allow_pickle=True, encoding="latin1"
    ).item()
video_object.update_paths(video_object_path)
## Load trajectories
trajectories = np.load(
    trajectories_path, allow_pickle=True, encoding="latin1"
).item()["trajectories"]

########################
# Set colours
########################

pal_dict_hex = {
    'iCab' : '#F1BB7B',
    'HdrR' : '#FA796C',
    'HNI' : '#AC3E3F',
    'Kaga' : '#79301F',
    'HO5' : '#D67236',
    'iCab_ref' : '#F1BB7B',
    'iCab_test' : '#AB7535'
    }

# Convert to RGB
def hex2rgb(hex):
    hex = hex.lstrip('#')
    rgb = tuple(int(hex[i:i+2], 16) for i in (0, 2, 4)) + (0,)
    return rgb

# Convert to RGB
def hex2bgr(hex):
    hex = hex.lstrip('#')
    bgr = tuple(int(hex[i:i+2], 16) for i in (4, 2, 0)) + (0,)
    return bgr

bgr_values = [hex2bgr(i) for i in pal_dict_hex.values()]
pal_dict_bgr = {list(pal_dict_hex.keys())[i] : bgr_values[i] for i in range(len(pal_dict_hex))}

########################
# Assign line ID and colours to trajectories object
########################

# Remove all frames with NA from the start
na_list = []
# Get list with "False" for all frames with NA
for i in range(len(trajectories)):
    na_list.append(np.isnan(trajectories[i]).any() == False)
## Filter trajectories file so that first frame is not NaN
first_frame = trajectories[na_list][0]

if REF_LOC == "Left" or "Right":
    # Pull out x coordinates
    x_coords = [first_frame[0][0], first_frame[1][0]] 
    # Get column index with min or max value based on REF location
    if REF_LOC == "Left":
        ref_col = np.argmin(x_coords)
    elif REF_LOC == "Right":
        ref_col = np.argmax(x_coords)

if REF_LOC == "Top" or "Bottom":
    # Pull out y coordinates
    y_coords = [first_frame[0][1], first_frame[1][1]] 
    # Get column index with min or max value based on REF location
    if REF_LOC == "Top":
        ref_col = np.argmin(y_coords)
    elif REF_LOC == "Bottom":
        ref_col = np.argmax(y_coords)

# If both fishes are from the same line, make the first object the reference 
if REF_LOC == "NA":
    ref_col = 0

# Create `cur_id_str` and `colors`` lists
if TEST_FISH == "iCab":
    cur_id_list = ["iCab_ref", "iCab_test"]
    colors = [pal_dict_bgr["iCab_ref"], pal_dict_bgr["iCab_test"]]
else:
    if ref_col == 0:
        cur_id_list = ["iCab", TEST_FISH]
        colors = [pal_dict_bgr["iCab"], pal_dict_bgr[TEST_FISH]]
    elif ref_col == 1:
        cur_id_list = [TEST_FISH, "iCab"]
        colors = [pal_dict_bgr[TEST_FISH], pal_dict_bgr["iCab"]]


# From `generate_trajectories_video()`

## args
centroid_trace_length=10
starting_frame=0
ending_frame=None

#ending_frame=300 # HASH OUT

## Get colours
#colors = get_spaced_colors_util(
#    video_object.number_of_animals, black=False
#)



## Set up fourcc
fourcc = cv2.VideoWriter_fourcc(*"XVID")

## Set up video writer
video_writer = cv2.VideoWriter(
    path_to_save_video,
    fourcc,
    video_object.frames_per_second,
    (video_object.width, video_object.height),
)

## Check whether there is a path to video segments (UNNECESSARY?)
if video_object.paths_to_video_segments is None:
    cap = cv2.VideoCapture(video_object.video_path)
else:
    cap = None

## Check whether there is an ending frame (UNNECESSARY?)
if ending_frame is None:
    ending_frame = len(trajectories)

## For each frame in the video
for frame_number in tqdm(
    range(starting_frame, ending_frame),
    desc="Generating video with trajectories...",
    ):

        #frame_number = 0 # HASH OUT
        # From `apply_func_on_frame()`

        ##
        segment_number = video_object.in_which_episode(frame_number)
        current_segment = segment_number

        if cap is None:      # (UNNECESSARY?)
            cap = cv2.VideoCapture(
                video_object.paths_to_video_segments[segment_number]
            )
            start = video_object._episodes_start_end[segment_number][0]
            cap.set(1, frame_number - start)
        else:
            cap.set(1, frame_number)
        ret, frame = cap.read()

        # Reduce video size 
        if video_object.resolution_reduction != 1:
            frame = cv2.resize(
                frame,
                None,
                fx=video_object.resolution_reduction,
                fy=video_object.resolution_reduction,
                interpolation=cv2.INTER_AREA,
            )        

        # Run `writeIds()` over each frame

        ordered_centroid = trajectories[frame_number]
        font = cv2.FONT_HERSHEY_SIMPLEX
        font_size = 0.5 * video_object.resolution_reduction
        font_width = int(1 * video_object.resolution_reduction)
        font_width = 1 if font_width == 0 else font_width
        circle_size = int(2 * video_object.resolution_reduction)

        for cur_id, centroid in enumerate(ordered_centroid):
            if sum(np.isnan(centroid)) == 0:
                if frame_number > centroid_trace_length:
                    centroids_trace = trajectories[
                        frame_number - centroid_trace_length : frame_number, cur_id
                    ]
                else:
                    centroids_trace = trajectories[:frame_number, cur_id]
                cur_id_str = cur_id_list[cur_id]
                int_centroid = np.asarray(centroid).astype("int")
                cv2.circle(
                    frame, tuple(int_centroid), circle_size, colors[cur_id], -1
                )
                cv2.putText(
                    frame,
                    cur_id_str,
                    tuple(int_centroid),
                    font,
                    font_size,
                    colors[cur_id],
                    font_width,
                )

                for centroid_trace in centroids_trace:
                    if sum(np.isnan(centroid_trace)) == 0:
                        int_centroid = np.asarray(centroid_trace).astype("int")
                        cv2.circle(
                            frame,
                            tuple(int_centroid),
                            circle_size,
                            colors[cur_id],
                            -1,
                        )

        video_writer.write(frame)

video_writer.release()
