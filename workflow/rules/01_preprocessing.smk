rule split_videos:
    input:
        os.path.join(config["raw_data_dir"], "{sample}.avi")
    output:
        os.path.join(config["lts_dir"], "split/{assay}/{sample}_{quadrant}.mp4")
    log:
        os.path.join(config["working_dir"], "logs/{sample}/{assay}/{quadrant}.log")
    params:
        samples_file = config["samples_file"]
    container:
        config["opencv"]
    script:
        "../scripts/split_videos.py"