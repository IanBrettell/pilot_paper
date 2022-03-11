# pull tracking parameters from config/samples.csv
def get_vid_length(wildcards):
    if wildcards.assay == "open_field":
        start = samples_df.loc[samples_df["sample"] == wildcards.sample, "of_start"]
        end = samples_df.loc[samples_df["sample"] == wildcards.sample, "of_end"]
    elif wildcards.assay == "novel_object":
        start = samples_df.loc[samples_df["sample"] == wildcards.sample, "no_start"]
        end = samples_df.loc[samples_df["sample"] == wildcards.sample, "no_end"]
    vid_length = int(end) - int(start)
    return(vid_length)

def get_bgsub(wildcards):
    if wildcards.assay == "open_field":
        target_col = "bgsub_" + "of_" + wildcards.quadrant
        bgsub = samples_df.loc[samples_df["sample"] == wildcards.sample, target_col].values[0]
    elif wildcards.assay == "novel_object":
        target_col = "bgsub_" + "no_" + wildcards.quadrant
        bgsub = samples_df.loc[samples_df["sample"] == wildcards.sample, target_col].values[0]
    return(bgsub)

def get_intensity_floor(wildcards):
    if wildcards.assay == "open_field":
        target_col = "intensity_floor_" + "of_" + wildcards.quadrant
        int_floor = int(samples_df.loc[samples_df["sample"] == wildcards.sample, target_col])
    elif wildcards.assay == "novel_object":
        target_col = "intensity_floor_" + "no_" + wildcards.quadrant
        int_floor = int(samples_df.loc[samples_df["sample"] == wildcards.sample, target_col])
    return(int_floor)

def get_intensity_ceiling(wildcards):
    if wildcards.assay == "open_field":
        target_col = "intensity_ceiling_" + "of_" + wildcards.quadrant
        int_ceiling = int(samples_df.loc[samples_df["sample"] == wildcards.sample, target_col])
    elif wildcards.assay == "novel_object":
        target_col = "intensity_ceiling_" + "no_" + wildcards.quadrant
        int_ceiling = int(samples_df.loc[samples_df["sample"] == wildcards.sample, target_col])
    return(int_ceiling)

def get_area_floor(wildcards):
    if wildcards.assay == "open_field":
        target_col = "area_floor_" + "of_" + wildcards.quadrant
        area_floor = int(samples_df.loc[samples_df["sample"] == wildcards.sample, target_col])
    elif wildcards.assay == "novel_object":
        target_col = "area_floor_" + "no_" + wildcards.quadrant
        area_floor = int(samples_df.loc[samples_df["sample"] == wildcards.sample, target_col])
    return(area_floor)

def get_area_ceiling(wildcards):
    if wildcards.assay == "open_field":
        target_col = "area_ceiling_" + "of_" + wildcards.quadrant
        area_ceiling = int(samples_df.loc[samples_df["sample"] == wildcards.sample, target_col])
    elif wildcards.assay == "novel_object":
        target_col = "area_ceiling_" + "no_" + wildcards.quadrant
        area_ceiling = int(samples_df.loc[samples_df["sample"] == wildcards.sample, target_col])
    return(area_ceiling)

# adapt memory usage for tracking videos
def get_mem_mb(wildcards, attempt):
    return attempt * 10000

# Track videos with idtrackerai
## Note: `output` is set as `trajectories.npy` instead of `trajectories_wo_gaps.npy`, presumably because
## in videos where there are no crossovers, the latter file is not produced.
rule track_videos:
    input:
        rules.split_videos.output
    output:
        os.path.join(
            config["working_dir"],
            "split/{assay}/session_{sample}_{quadrant}/trajectories/trajectories.npy"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/track_videos/{assay}/{sample}/{quadrant}.log"
        ),
    params:
        vid_length = get_vid_length,
        vid_name = "{sample}_{quadrant}",
        bgsub = get_bgsub,
        intensity_floor = get_intensity_floor,
        intensity_ceiling = get_intensity_ceiling,
        area_floor = get_area_floor,
        area_ceiling = get_area_ceiling,
    resources:
        mem_mb = 10000
    container:
        config["idtrackerai"]
    shell:
        """
        idtrackerai terminal_mode \
            --_video {input} \
            --_bgsub '{params.bgsub}' \
            --_range [0,{params.vid_length}] \
            --_nblobs 2 \
            --_intensity [{params.intensity_floor},{params.intensity_ceiling}] \
            --_area [{params.area_floor},{params.area_ceiling}] \
            --_session {params.vid_name} \
            --exec track_video \
                2> {log}
        """

def get_trajectories_file(wildcards):
    # Get path of trajectories files
    traj_wo_gaps_file = os.path.join(
        config["working_dir"],
        "split/{assay}/session_{sample}_{quadrant}/trajectories_wo_gaps/trajectories_wo_gaps.npy")
    traj_file = os.path.join(
        config["working_dir"],
        "split/{assay}/session_{sample}_{quadrant}/trajectories/trajectories.npy")
    # If there is no `trajectories_wo_gaps.npy` file, return the `trajectories.npy` file
    if os.path.exists(traj_wo_gaps_file):
        return(traj_wo_gaps_file)
    else:
        return(traj_file)

# Convert numpy arrays to .csv files
# Not needed if `CONVERT_TRAJECTORIES_DICT_TO_CSV_AND_JSON = True`
# is added to `local_settings.py`
#rule trajectories_to_csv:
#    input:
#        trajectories = rules.track_videos.output,
#        script = "workflow/scripts/trajectories_to_csv.py"
#    output:
#        os.path.join(config["working_dir"], "split/{assay}/session_{sample}_{quadrant}/trajectories/trajectories.trajectories.csv")
#    log:
#        os.path.join(config["working_dir"], "logs/trajectories_to_csv/{assay}/{sample}/{quadrant}.log"),
#    params:
#        in_path = os.path.join(config["working_dir"], "split/{assay}/session_{sample}_{quadrant}")
#    resources:
#        mem_mb = 100,
#    shell:
#        """
#        python {input.script} {params.in_path} \
#            2> {log}
#        """