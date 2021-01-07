#!/bin/sh

# Get variables from script

#input=$(echo $1)
#intensity=$(echo $2)
#area=$(echo $3)
#range=$(echo $4)
#session_name=$(echo $5)

samples_file=$(echo $1)
input=$(echo $2)
assay=$(basename $input | cut -f7-8 -d"_" | sed 's/.mp4//')
run=$(basename $input | cut -f1-5 -d"_")
target_line=$(grep $run $samples_file)
intensity=$(echo "[$(echo $target_line | cut -f16 -d','),$(echo $target_line | cut -f17 -d',')]")
area=$(echo "[$(echo $target_line | cut -f18 -d','),$(echo $target_line | cut -f19 -d',')]")
if [ "$assay" = "novel_object" ]; then
  range=$(echo "[0,$(echo $target_line | cut -f15 -d',')]")
else [ "$assay" = "open_field" ]
  range=$(echo "[0,$(echo $target_line | cut -f14 -d',')]")
fi
session_name=$(basename $input | cut -f1 -d".")

## Activate conda env
#
#source activate idtrackerai

# Run idtrackerai

idtrackerai terminal_mode \
  --_video $input \
  --_bgsub 'True' \
  --_intensity $intensity \
  --_area $area \
  --_range $range \
  --_nblobs 2 \
  --_session $session_name \
  --exec track_video
