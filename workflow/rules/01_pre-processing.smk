## Copy videos from FTP to Codon (for ease downstream)
#rule copy_videos:
#    input:
#        os.path.join(config["raw_data_dir"], "{sample}.avi"),
#    output:
#        os.path.join(config["working_dir"], "raw_videos/{sample}.avi"),
#    log:
#        os.path.join(config["working_dir"], "logs/copy_videos/{sample}.log"),
#    shell:
#        """
#        cp {input} {output} \
#            2> {log}
#        """

# Generate single-frame grab showing coordinate of splits
rule set_split_coords:
    input:
        video = os.path.join(config["raw_data_dir"], "{sample}.avi"),
    output:
        fig = "results/split_coord_images/{sample}.png",
    log:
        os.path.join(config["working_dir"], "logs/set_split_coords/{sample}.log")
    params:
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
        os.path.join(config["data_store_dir"], "split/{assay}/{sample}_{quadrant}.mp4"),
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


