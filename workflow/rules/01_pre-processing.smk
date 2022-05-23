# When tank_side == "R", video is flipped 180 degrees (camera was upside-down)
# This rule flips those videos back in the process of creating .avi copies of all videos
rule adjust_orientation:
    input:
        os.path.join(config["raw_data_dir"], "{sample}.avi"),
    output:
        os.path.join(config["working_dir"], "flipped/{sample}.avi"),
    log:
        os.path.join(config["working_dir"], "logs/adjust_orientation/{sample}.log")
    params:
        sample = "{sample}",
        samples_file = lambda wildcards: config["samples_file"]
    container:
        config["opencv"]
    resources:
        mem_mb = 1000
    script:
        "../scripts/adjust_orientation.py"

# Generate single-frame grab showing coordinate of splits
rule set_split_coords:
    input:
        video = rules.adjust_orientation.output,
    output:
        fig = "results/split_coord_images/{assay}/{sample}.png",
    log:
        os.path.join(config["working_dir"], "logs/set_split_coords/{assay}/{sample}.log")
    params:
        assay = "{assay}",
        sample = "{sample}",
        samples_file = lambda wildcards: config["samples_file"]
    container:
        config["opencv"]
    resources:
        mem_mb = 500
    script:
        "../scripts/set_split_coords.py"

# Split videos into quadrants and assays (1 raw video * 4 quadrants * 2 assays = 8 output videos)
rule split_videos:
    input:
        rules.adjust_orientation.output,
    output:
        os.path.join(
            config["working_dir"],
            "split/{assay}/{sample}_{quadrant}.avi"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/split_videos/{assay}/{sample}/{quadrant}.log"
        ),
    params:
        sample = "{sample}",
        assay = "{assay}",
        quadrant = "{quadrant}",
        samples_file = lambda wildcards: config["samples_file"]
    container:
        config["opencv"]
    resources:
        mem_mb = 1000
    script:
        "../scripts/split_videos.py"

rule get_split_video_dims:
    input:
        expand(rules.split_videos.output,
                zip,
                assay = ASSAYS,
                sample = SAMPLES,
                quadrant = QUADRANTS        
        ),
    output:
        "config/split_video_dims.csv"
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_split_video_dims/all.log"
        ),
    container:
        config["opencv"]
    resources:
        mem_mb = 300
    script:
        "../scripts/get_split_video_dims.py"            


