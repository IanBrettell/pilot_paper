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
    fps = int(samples_df.loc[samples_df["sample"] == wildcards.sample, "fps"])
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
            config["data_store_dir"],
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

# Generate videos with coloured trails superimposed
rule coloured_trails:
    input:
        video_object=os.path.join(
            config["working_dir"],
            "split/{assay}/session_{sample}_{quadrant}/video_object.npy",
        ),
        trajectories=get_trajectories_file,
    output:
        os.path.join(
            config["working_dir"],
            "split/{assay}/{sample}_{quadrant}_tracked.avi",
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/coloured_trails/{assay}/{sample}/{quadrant}.log"
        ),
    container:
        config["idtrackerai"]
    resources:
        mem_mb=5000,
    script:
        "../scripts/coloured_trails.py"
    
