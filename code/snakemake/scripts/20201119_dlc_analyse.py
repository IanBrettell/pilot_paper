#!/usr/bin/env python3

# Import libaries

import os
os.environ["DLClight"]="True" # suppress GUI

import deeplabcut

# Setup project variables

project_path = snakemake.params.dlc_project_path
#config_path = os.path.join(snakemake.params.dlc_project_path, 'config.yaml')
# dlc_project_path = 'dlc/pp_proj-ian_brettell-2020-11-17'
# config_path = os.path.join(dlc_project_path, 'config.yaml')

target_video = snakemake.input[0]
# target_video = 'dlc/pp_proj-ian_brettell-2020-11-17/videos/20190611_1410_icab_kaga_R_q4_open_field.mp4'
# get basename for name of copied folder
target_video_bname = os.path.basename(target_video).strip('.mp4')

# Make functions to copy project folder
# (with assistance from: https://www.pythoncentral.io/how-to-recursively-copy-a-directory-folder-in-python/)

import shutil
import errno

to_ignore = 'videos'

def ignore_function(ignore):
    def _ignore_(path, names):
        ignored_names = []
        if ignore in names:
            ignored_names.append(ignore)
        return set(ignored_names)
    return _ignore_

def copy(src, dest):
    try:
        shutil.copytree(src, dest, ignore=ignore_function(to_ignore))
    except OSError as e:
        # If the error was caused because the source wasn't a directory
        if e.errno == errno.ENOTDIR:
            shutil.copy(src, dest)
        else:
            print('Directory not copied. Error: %s' % e)

# Copy dlc project folder minus the videos folder

tmp_dir = os.path.join('tmp', target_video_bname) # create path of tmp directory
copy(src = project_path, dest = tmp_dir) # copy
os.mkdir(os.path.join(tmp_dir, 'videos')) # make video dir

# Change project path in tmp config file

config_path = os.path.join(tmp_dir, 'config.yaml')

new_line = 'project_path: ' + tmp_dir + '\n'# create new line

with open(config_path, 'r') as file:
    # read a list of lines into data
    data = file.readlines()

data[7] = new_line # substitute 8th line

with open(config_path, 'w') as file: # write to new file
    file.writelines( data )

# Analyse video

deeplabcut.analyze_videos(config_path, [target_video], TFGPUinference = False)

# Copy output to original folder

from distutils.dir_util import copy_tree
copy_tree(os.path.join(tmp_dir, 'videos'), os.path.join(project_path, 'videos'))

# Delete tmp project folder

shutil.rmtree(tmp_dir)
