# Import tensorflow, suppress GUI, import DeepLabCut

import sys

import tensorflow
print(tensorflow.__version__)

import os
os.environ["DLClight"]="True"

import deeplabcut
print(deeplabcut.__version__)

# Get inputs

config_path = sys.argv[1]
#config_path = 'dlc/pp_proj-ian_brettell-2020-12-03/config.yaml'
target_video = sys.argv[2]
#videofile_path = ['dlc/pp_proj-ian_brettell-2020-12-03/videos/20190613_1617_icab_ho5_L_q1_open_field.mp4']

# Analyse videos

deeplabcut.analyze_videos(config_path, videofile_path, videotype='mp4')
