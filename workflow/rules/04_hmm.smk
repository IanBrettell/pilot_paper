rule movement_metrics:
    input:
        rules.assign_ref_test.output,
    output:
        os.path.join(
            config["working_dir"],
            "with_metrics/{assay}/{sample}/{quadrant}/{interval}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/movement_metrics/{assay}/{sample}/{quadrant}/{interval}.log"
        ),
    params:
        seconds_interval = "{interval}",
        source_file = "workflow/scripts/movement_metrics_source.R",
    resources:
        mem_mb = 500
    container:
        config["rocker_tidyverse"]
    script:
        "../scripts/movement_metrics.R"

def get_hmm_variables(wildcards):
    

rule hmm:
    input:
        rules.movement_metrics.output,
    output:
        os.path.join(
            config["working_dir"],
            "hmm_out/{assay}/{sample}/{quadrant}/{interval}/{variables}/{n_states}.csv"
        ),        
    log:
        os.path.join(
            config["working_dir"],
            "logs/hmm/{assay}/{sample}/{quadrant}/{interval}/{variables}/{n_states}.log"
        ),
    params:
        n_states = "{n_states}",
        variables = "{variables}"
    resources:
        mem_mb = 1000,
    container:
        config["hmmlearn"]
    script:
        "../scripts/hmm.py"
