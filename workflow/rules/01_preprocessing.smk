rule copy_videos:
    input:
        os.path.join(config["raw_data_dir"], "{sample}.avi"),
    output:
        os.path.join(config["working_dir"], "raw_videos/{sample}.avi"),
    log:
        os.path.join(config["working_dir"], "logs/copy_videos/{sample}.log"),
    shell:
        """
        cp {input} {output} \
            2> {log}
        """

rule split:
    input:
        os.path.join(config["input_dir"], "{sample}.avi")
    params:
        samples_file = config["samples_file"]
    output:
        "split/{assay}/{sample}_{quadrant}_{assay}.mp4"
    conda:
        config["python_env"]
    script:
        config["split_script"]