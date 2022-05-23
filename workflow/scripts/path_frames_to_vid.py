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
from plotnine import *
import os
import shutil

# Get variables

## Debug
#IN = "/hps/nobackup/birney/users/ian/pilot/hmm_out/0.08/dist_angle/15.csv"
#DIMS = "config/split_video_dims.csv"
#ASSAY = "novel_object"
#SAMPLE = "20190613_0953_icab_hni_R"
#FPS = 30
#TMP = "/hps/nobackup/birney/users/ian/pilot/tmp"

## True
IN = snakemake.input.hmm[0]
DIMS = snakemake.input.dims[0]
OUT = snakemake.output[0]
ASSAY = snakemake.params.assay
SAMPLE = snakemake.params.sample
FPS = int(snakemake.params.fps)
TMP = snakemake.resources.tmpdir

#######################
# Set plotting parameters
#######################

REF = SAMPLE.split('_')[2]
TEST = SAMPLE.split('_')[3]

if REF == TEST:
    line_dict = {
        "ref" : "iCab\nref",
        "test" : "iCab\ntest"
    }
    pal_dict = {
        "iCab\nref" : "#F1BB7B",
        "iCab\ntest" : "#AB7535"
    }
else:
    line_dict = {
        "icab" : "iCab",
        "hdr" : "HdrR",
        "hni" : "HNI",
        "kaga" : "Kaga",
        "ho5": "HO5"
    }
    pal_dict = {
        "iCab" : "#F1BB7B",
        "HdrR" : "#FA796C",
        "HNI" : "#AC3E3F",
        "Kaga" : "#79301F",
        "HO5" : "#D67236"
    }

#######################
# Read in tracking data
#######################

# Get date and time
DATE = SAMPLE.split('_')[0]
TIME = SAMPLE.split('_')[1]

# Read in file

df = pd.read_csv(IN)
# Convert time to string
df['time'] = df['time'].astype('string')
# Add a 0 to the start of `time` if only 3 characters
def add_0(x):
    if len(x) == 4:
        return x
    elif len(x) == 3:
        return "0" + x

df['time'] = df['time'].apply(add_0)

# Filter for assay and sample
df_filt = df.loc[(df['assay'] == ASSAY) & (df['date'] == int(DATE)) & (df['time'] == TIME)].copy()

# Recode states

old_states = df.groupby(['state'])['distance'].mean().sort_values().index
new_states = list(range(1, len(old_states) + 1))

state_dict = {}
for x,y in zip(old_states, new_states):
    state_dict[x] = y

df_filt["state_recode"] = df_filt["state"].map(state_dict)

# Create `line` column and recode

if REF == TEST:
    df_filt['line'] = df_filt['fish'].map(line_dict)
else:
    conditions = [
        (df_filt['fish'] == "ref"),
        (df_filt['fish'] == "test")
    ]
    values = [
        df_filt['ref_fish'],
        df_filt['test_fish']
    ]
    df_filt['line'] = np.select(conditions, values)
    df_filt['line'] = df_filt['line'].map(line_dict)

#######################
# Read in dims
#######################

dims = pd.read_csv(DIMS)
dims = dims.loc[(dims['sample'] == SAMPLE) & (dims['assay'] == ASSAY)]
N_FRAMES = max(dims['n_frames'])

# Set order of quadrants
df_filt['quadrant'] = pd.Categorical(df_filt['quadrant'], categories = ["q2", "q1", "q3", "q4"])

# Get max width and height

wid = dims['wid'].max()
hei = dims['hei'].max()

# Get total height and width in pixels (to match with the raw videos)

TOT_WID = dims.loc[dims['quadrant'].isin(['q1', 'q2'])]['wid'].sum() 
TOT_HEI = dims.loc[dims['quadrant'].isin(['q1', 'q4'])]['hei'].sum() 

####################
# Set up video writer
####################

fourcc = cv2.VideoWriter_fourcc('h', '2', '6', '4')

## Debug
#video_writer = cv2.VideoWriter(
#    "/hps/nobackup/birney/users/ian/pilot/tmp_out.avi",
#    fourcc,
#    FPS,
#    (TOT_WID, TOT_HEI),
#    isColor = True
#)

video_writer = cv2.VideoWriter(
    OUT,
    fourcc,
    FPS,
    (TOT_WID, TOT_HEI),
    isColor = True
)

#######################
# Plot frames and write to video
#######################

# Get all available frames
all_frames = list(range(1, N_FRAMES + 1))

# Get all frames we have
rec_frames = df_filt['frame'].unique().tolist()
rec_frames.sort()

for i in all_frames[0:500]:
    # If the frame is not included in the frames we have data for...
    if i not in rec_frames:
        # get the next frame we do have
        filt_frames = []
        for frame in rec_frames:
            if frame > i:
                filt_frames.append(frame)
        plot_frame = min(filt_frames)
    elif i in rec_frames:
        plot_frame = i
    print(plot_frame)
    # If file already exists, write directly to file
    out_path = os.path.join(
        TMP,
        "path_frames",
        ASSAY,
        SAMPLE, 
        str(plot_frame) + ".png"
    )
    # Make directory
    os.makedirs(os.path.dirname(out_path), exist_ok = True)
    # If the plot .png is already there, read it in and write
    if os.path.exists(out_path):
        # read plot
        frame = cv2.imread(out_path)
        # resize
        frame_out = cv2.resize(frame, (TOT_WID, TOT_HEI))
        # write
        video_writer.write(frame_out)
    # otherwise create the plot
    else:
        ## Filter df
        dat = df_filt.loc[df_filt['frame'] <= plot_frame]
        ## Plot
        plot = (ggplot(dat) +
                geom_path(aes('x', 'y', colour = 'line')) +
                facet_wrap('quadrant', nrow = 2) +
                scale_color_manual(pal_dict) +
                scale_x_continuous(limits = [0,wid]) +
                scale_y_reverse(limits = [hei,0]) +
                guides(color = None) +
                theme(
                    aspect_ratio = 1,
                    strip_background = element_blank(),
                    strip_text = element_blank(),
                    axis_line = element_blank(),
                    axis_text = element_blank(),
                    axis_title = element_blank(),
                    axis_ticks = element_blank(),
                    panel_background = element_blank(),
                    plot_background = element_rect(fill = "white")
                )
        )
        # Save
        #out_path = os.path.join(TMP, SAMPLE, REF_TEST, str(plot_frame) + ".png")
        plot.save(out_path, width = 9, height = 9)
        # Write to video
        ## read image
        frame = cv2.imread(out_path)
        # resize
        frame_out = cv2.resize(frame, (TOT_WID, TOT_HEI))
        # write
        video_writer.write(frame_out)
        # delete file
        os.remove(out_path)

video_writer.release()

# remove temp folder

shutil.rmtree(os.path.dirname(out_path))
