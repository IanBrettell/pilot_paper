
conda activate idtrackerai

## Sample parameters file
samples_file=pilot_paper/data/20210118_paramaters.csv
## Video store
video_store=gs://medaka-video-store

for in_sample in $(gsutil ls $video_store/videos/20190613* | cut -f5 -d'/' | cut -f1 -d'.' ); do
  # check whether final trajectories_wo_gaps.npy file exists in results folder
  gsutil -q stat $video_store/results/session_$in_sample/trajectories/trajectories.npy ;
  status=$? ;
  # Run `idtrackerai` if results don't exist
  if [[ $status == 1 ]]; then
    # Create log
    echo "$in_sample: started $(date)" >> pilot_paper/data/20210211_test_train.txt ;
    # Set training folder
#    target_line=$(echo $in_sample | cut -f4 -d'_' )
#    echo "KNOWLEDGE_TRANSFER_FOLDER_IDCNN = trained/$target_line" > local_settings.py ;
    echo "KNOWLEDGE_TRANSFER_FOLDER_IDCNN = accumulation_0" > local_settings.py ;
    # Copy video to VM
    gsutil -m cp $video_store/videos/$in_sample.mp4 videos/ ;
    # Get parameters
    input_video=$(echo videos/$in_sample.mp4 ) ;
    target_string=$( grep $in_sample $samples_file ) ;
    vid_length=$( echo $target_string | cut -f2 -d',' ) ;
    int_floor=$( echo $target_string | cut -f3 -d',' ) ;
    int_ceiling=$( echo $target_string | cut -f4 -d',' ) ;
    area_floor=$( echo $target_string | cut -f5 -d',' ) ;
    area_ceiling=$( echo $target_string | cut -f6 -d',' ) ;
    # Run `idtrackerai`
    idtrackerai terminal_mode \
              --_video $input_video \
              --_bgsub 'True' \
              --_intensity [$int_floor,$int_ceiling] \
              --_area [$area_floor,$area_ceiling] \
              --_range [0,$vid_length] \
              --_nblobs 2 \
              --_session $in_sample \
              --exec track_video ;
    # Run script to convert tracking data to CSVs
    python pilot_paper/code/scripts/trajectories_to_csv.py videos/session_$in_sample ;
    # Copy to storage bucket
    gsutil -m cp -r \
      videos/session_$in_sample \
      gs://medaka-video-store/results ;
    # Remove video and results
    rm videos/$in_sample.mp4
    rm -rf videos/session_$in_sample
    # Confirm with output
    echo "$in_sample: completed $(date)" >> pilot_paper/data/20210211_test_train.txt
    echo "$in_sample: VIDEO PROCESSED AND COPIED TO BUCKET. PROCEEDING TO NEXT..."
  else
    echo "$in_sample: VIDEO ALREADY PROCESSED. PROCEEDING TO NEXT..."
  fi
done
