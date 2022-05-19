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

# Get variables

IN = "/hps/nobackup/birney/users/ian/pilot/hmm_out/0.08/dist_angle/15.csv"
ASSAY = "novel_object"
SAMPLE = "20190613_0953_icab_hni_R"
REF_TEST = "test"
DIMS = "config/split_video_dims.csv"
TMP = "/hps/nobackup/birney/users/ian/pilot/tmp"
FPS = 30



#######################
# Set plotting parameters
#######################

if REF_TEST == "test":
    pal_option = "viridis"
elif REF_TEST == "ref":
  pal_option = "inferno"

#######################
# Read in tracking data
#######################

# Create line recode vector
line_dict = {
    "icab" : "iCab",
    "hdr" : "HdrR",
    "hni" : "HNI",
    "kaga" : "Kaga",
    "ho5": "HO5"
}

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

# Filter for sample
df_filt = df.loc[(df['date'] == int(DATE)) & (df['time'] == TIME)]

# Recode states

old_states = df.groupby(['state'])['distance'].mean().sort_values().index
new_states = list(range(1, len(old_states) + 1))

state_dict = {}
for x,y in zip(old_states, new_states):
    state_dict[x] = y

df_filt["state_recode"] = df_filt["state"].map(state_dict)

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
video_writer = cv2.VideoWriter(
    "/hps/nobackup/birney/users/ian/pilot/tmp_out.avi",
    fourcc,
    FPS,
    (TOT_WID, TOT_HEI),
    isColor = True
)

#######################
# Create a dictionary of missing frames
#######################

# Get all available frames
all_frames = list(range(1, N_FRAMES + 1))

# Get all frames we have
rec_frames = df_filt['frame'].unique().tolist()
rec_frames.sort()

for i in all_frames[0:10]:
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
    out_path = os.path.join(TMP, SAMPLE, REF_TEST, plot_frame + ".png")
    if os.path.exists(out_path):
        frame = cv2.imread(out_path)
        video_writer.write(frame)
    # otherwise create the plot
    else:
        ## Filter df
        dat = df_filt.loc[df_filt['frame'] <= plot_frame]
        ## Plot
        plot = (ggplot(dat) +
                geom_point(aes('x', 'y', colour = 'state_recode'),
                            alpha = 0.8) +
                facet_wrap('quadrant', nrow = 2) +
                scale_color_continuous(pal_option) +
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
        out_path = os.path.join(TMP, SAMPLE, REF_TEST, plot_frame + ".png")
        plot.save(out_path, width = 9, height = 9)
        # Write to video
        ## read image
        frame = cv2.imread(out_path)
        video_writer.write(frame)

video_writer.release()











all_frames = pd.DataFrame({"frame" : list(range(1, N_FRAMES + 1))})
all_frames = all_frames.merge(df_filt.groupby(['quadrant']), on = 'frame', how = 'left')
# Fill forward
all_frames = all_frames.fillna(method = "bfill")

QUADS = ["q2", "q1", "q3", "q4"]
df_dict = {}
for quad in QUADS:
    quad_df = df_filt.loc[df_filt['quadrant'] == quad]
    all_df = all_frames.merge(quad_df, on = 'frame', how = 'left')
    all_df.fillna(method = "bfill")
    df_dict[quad] = all_df

#######################
# Plot
#######################



#for FRAME in range(0, 1):
FRAME = 0
dat = df_filt.loc[df_filt['frame'] <= 8000]
plot = (ggplot(dat) +
        geom_point(aes('x', 'y', colour = 'state_recode'),
                    alpha = 0.8) +
        facet_wrap('quadrant', nrow = 2) +
        scale_color_continuous(pal_option) +
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
plot.save("p9.png", width = 9, height = 9)
