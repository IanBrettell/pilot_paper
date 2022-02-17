# When tank_side == "R", video is flipped 180 degrees (camera was upside-down)
# This rule flips those videos back in the process of creating .mp4 copies of all videos
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
    script:
        "../scripts/adjust_orientation.py"

# Decompress videos so that they can be viewed in Fiji
# in order to determine the start and end frames for open_field and novel_object assays
# NOTE: creates huge files, so consider hashing out the rule and deleting the files once
# metadata has been collected
rule decompress_videos:
    input:
        rules.adjust_orientation.output,
    output:
        os.path.join(config["working_dir"], "decompressed/{sample}.avi"),
    log:
        os.path.join(config["working_dir"], "logs/decompress_videos/{sample}.log"),
    container:
        config["ffmpeg"]
    shell:
        """
        ffmpeg -i {input} -an -vcodec rawvideo -y {output} \
            2> {log}
        """

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
    script:
        "../scripts/set_split_coords.py"

# Split videos into quadrants and assays (1 raw video * 4 quadrants * 2 assays = 8 output videos)
rule split_videos:
    input:
        rules.set_split_coords.input,
    output:
        os.path.join(config["data_store_dir"], "split/{assay}/{sample}_{quadrant}.avi"),
    log:
        os.path.join(config["working_dir"], "logs/split_videos/{assay}/{sample}/{quadrant}.log"),
    params:
        sample = "{sample}",
        assay = "{assay}",
        quadrant = "{quadrant}",
        samples_file = lambda wildcards: config["samples_file"]
    container:
        config["opencv"]
    script:
        "../scripts/split_videos.py"


