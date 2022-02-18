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


def writeIds(
    video_object,
    frame,
    frame_number,
    trajectories,
    centroid_trace_length,
    colors,
):
    ordered_centroid = trajectories[frame_number]
    font = cv2.FONT_HERSHEY_SIMPLEX
    font_size = 1 * video_object.resolution_reduction
    font_width = int(3 * video_object.resolution_reduction)
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
            cur_id_str = str(cur_id + 1)
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

    return frame


def apply_func_on_frame(
    video_object,
    trajectories,
    frame_number,
    video_writer,
    colors,
    cap=None,
    func=None,
    centroid_trace_length=10,
):
    segment_number = video_object.in_which_episode(frame_number)
    current_segment = segment_number

    if cap is None:
        cap = cv2.VideoCapture(
            video_object.paths_to_video_segments[segment_number]
        )
        start = video_object._episodes_start_end[segment_number][0]
        cap.set(1, frame_number - start)
    else:
        cap.set(1, frame_number)
    ret, frame = cap.read()

    if ret:
        if video_object.resolution_reduction != 1:
            frame = cv2.resize(
                frame,
                None,
                fx=video_object.resolution_reduction,
                fy=video_object.resolution_reduction,
                interpolation=cv2.INTER_AREA,
            )
        frame = func(
            video_object,
            frame,
            frame_number,
            trajectories,
            centroid_trace_length,
            colors,
        )

        return frame


def generate_trajectories_video(
    video_object,
    trajectories,
    path_to_save_video,
    func=writeIds,
    centroid_trace_length=10,
    starting_frame=0,
    ending_frame=None,
):

    colors = get_spaced_colors_util(
        video_object.number_of_animals, black=False
    )

    fourcc = cv2.VideoWriter_fourcc(*"XVID")

    video_writer = cv2.VideoWriter(
        path_to_save_video,
        fourcc,
        video_object.frames_per_second,
        (video_object.width, video_object.height),
    )

    if video_object.paths_to_video_segments is None:
        cap = cv2.VideoCapture(video_object.video_path)
    else:
        cap = None

    if ending_frame is None:
        ending_frame = len(trajectories)

    for frame_number in tqdm(
        range(starting_frame, ending_frame),
        desc="Generating video with trajectories...",
    ):
        frame = apply_func_on_frame(
            video_object,
            trajectories,
            frame_number,
            video_writer,
            colors,
            cap=cap,
            func=writeIds,
            centroid_trace_length=centroid_trace_length,
        )
        video_writer.write(frame)


def main(video_object_path, trajectories_path, path_to_save_video):
    video_object = np.load(
        video_object_path, allow_pickle=True, encoding="latin1"
    ).item()
    video_object.update_paths(video_object_path)
    trajectories = np.load(
        trajectories_path, allow_pickle=True, encoding="latin1"
    ).item()["trajectories"]

    generate_trajectories_video(
        video_object,
        trajectories,
        path_to_save_video,
    )


# NOTE: it is fundamental to convert snakemake objects to strings.
# opencv drops any input that is in a type that it doesn't expect.
# This causes very confusing errors that are hard to reconduce to a type error:
# the order of arguments in a function is altered and all sorts of things can
# happen.

main(
    video_object_path=str(snakemake.input.video_object),
    trajectories_path=str(snakemake.input.trajectories),
    path_to_save_video=str(snakemake.output),
)