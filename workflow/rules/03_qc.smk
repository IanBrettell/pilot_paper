# Choose input trajectories .csv file
def get_final_csvs(wildcards):
    # Get path of csv files
    ## Trajectories without gaps
    traj_wo_gaps_file = os.path.join(
        config["working_dir"],
        "split",
        wildcards.assay,
        "session_" + wildcards.sample + "_" + wildcards.quadrant,
        "trajectories_wo_gaps",
        "trajectories_wo_gaps.trajectories.csv"
        )
    ## Trajectories (with gaps)
    traj_file = os.path.join(
        config["working_dir"],
        "split",
        wildcards.assay,
        "session_" + wildcards.sample + "_" + wildcards.quadrant,
        "trajectories",
        "trajectories.trajectories.csv"
        )
    # If there is no "without gaps" file, return the "trajectories" file
    if os.path.exists(traj_wo_gaps_file):
        return(traj_wo_gaps_file)
    elif os.path.exists(traj_file):
        return(traj_file)

# Get frames-per-second
def get_fps(wildcards):
    fps = samples_df.loc[samples_df["sample"] == wildcards.sample, "fps"]
    return(fps)

# Get relative location of reference iCab fish
def get_ref_loc(wildcards):
    if wildcards.assay == "open_field":
        target_col = "cab_coords_" + "of_" + wildcards.quadrant
        ref_loc = samples_df.loc[samples_df["sample"] == wildcards.sample, target_col]
    elif wildcards.assay == "novel_object":
        target_col = "cab_coords_" + "no_" + wildcards.quadrant
        ref_loc = samples_df.loc[samples_df["sample"] == wildcards.sample, target_col]
    ref_loc = ref_loc.values[0]
    # `ref_loc` is nan if there are two iCabs, so convert to string
    if pd.isna(ref_loc):
        ref_loc = "NA"
    return(ref_loc)

# Assign reference and test fish IDs, and filter for frames up to 10 minutes
rule assign_ref_test:
    input:
        get_final_csvs,
    output:
        os.path.join(
            config["working_dir"],
            "final_tracks/{assay}/{sample}_{quadrant}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/assign_ref_test/{assay}/{sample}/{quadrant}.log"
        ),
    params:
        fps = get_fps,
        ref_loc = get_ref_loc,
    resources:
        mem_mb = 200
    script:
        "../scripts/assign_ref_test.py"
    
# Create .csv with proportion of frames tracked
rule tracking_success:
    input:
        expand(rules.assign_ref_test.output,
            zip,
            assay = ASSAYS,
            sample = SAMPLES,
            quadrant = QUADRANTS                      
        ),
    output:
        "config/tracking_success.csv"
    log:
        os.path.join(
            config["working_dir"],
            "logs/tracking_success/tracking_success.log"
        ),
    container:
        config["rocker_tidyverse"]
    resources:
        mem_mb = 1000
    script:
        "../scripts/tracking_success.R"

# Function to pull trajectories.npy if trajectories_wo_gaps.npy does not exist
def get_trajectories_file(wildcards):
    # Get path of csv files
    ## Trajectories without gaps
    traj_wo_gaps_file = os.path.join(
        config["working_dir"],
        "split",
        wildcards.assay,
        "session_" + wildcards.sample + "_" + wildcards.quadrant,
        "trajectories_wo_gaps",
        "trajectories_wo_gaps.npy"
        )
    ## Trajectories (with gaps)
    traj_file = os.path.join(
        config["working_dir"],
        "split",
        wildcards.assay,
        "session_" + wildcards.sample + "_" + wildcards.quadrant,
        "trajectories",
        "trajectories.npy"
        )
    # If there is no "without gaps" file, return the "trajectories" file
    if os.path.exists(traj_wo_gaps_file):
        return(traj_wo_gaps_file)
    elif os.path.exists(traj_file):
        return(traj_file)
    else:
        raise Exception('No trajectories file found.')

# Generate videos with coloured trails superimposed
#rule coloured_trails:
#    input:
#        video_object=os.path.join(
#            config["working_dir"],
#            "split/{assay}/session_{sample}_{quadrant}/video_object.npy",
#        ),
#        trajectories=get_trajectories_file,
#    output:
#        os.path.join(
#            config["working_dir"],
#            "tracked/{assay}/{sample}_{quadrant}.avi",
#        ),
#    log:
#        os.path.join(
#            config["working_dir"],
#            "logs/coloured_trails/{assay}/{sample}/{quadrant}.log"
#        ),
#    container:
#        config["idtrackerai"]
#    resources:
#        mem_mb=5000,
#    script:
#        "../scripts/coloured_trails.py"
    
rule coloured_trails_labels:
    input:
        video_object=rules.track_videos.output.video_obj,
        trajectories=get_trajectories_file,
    output:
        os.path.join(
            config["working_dir"],
            "tracked/{assay}/{sample}_{quadrant}.avi",
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/coloured_trails_labels/{assay}/{sample}/{quadrant}.log"
        ),
    params:
        sample = "{sample}",
        ref_loc = get_ref_loc,
    container:
        config["idtrackerai"]
    resources:
        mem_mb=5000,
    script:
        "../scripts/coloured_trails_labels.py"

rule stitch_tracked_vids:
    input:
        expand(os.path.join(
            config["working_dir"],
            "tracked/{{assay}}/{{sample}}_{quadrant}.avi"
            ),
                quadrant = QUADRANTS_ALL
        ),
    output:
        os.path.join(
            config["working_dir"],
            "stitched/{assay}/{sample}.avi",
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/stitch_tracked_vids/{assay}/{sample}.log"
        ),
    params:
        fps = get_fps,
    container:
        config["opencv"]
    resources:
        mem_mb=5000,
    script:
        "../scripts/stitch_tracked_vids.py"

rule pull_visual_check_sample:
    input:
        config["samples_file"],
    output:
        csv = "config/random_visual_check_sample.csv",
    log:
        os.path.join(
            config["working_dir"],
            "logs/pull_visual_check_sample/all.log"
        ),
    params:
        n_samples = 20,
        seed = 5,
        seed_assay = 10
    resources:
        mem_mb = 200,
    run:
        # Set up log
        with open(log[0], "w") as f:
            sys.stderr = sys.stdout = f
        import random
        # Read in samples file
        df = pd.read_csv(input[0])
        # Pull out sample names
        VIDEOS = df['sample'].values
        # Get random sample of videos
        random.seed(params.seed)
        # NOTE: need to sort VIDEOS because `random` only takes a sequence
        rand_vids = sorted(random.sample(sorted(VIDEOS), params.n_samples))
        # Create DF
        out = pd.DataFrame({'sample': rand_vids})
        # Write to file
        out.to_csv(output[0], index=False)
