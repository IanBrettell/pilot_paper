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

rule merge_csvs:
    input:
        expand(os.path.join(
            config["working_dir"],
            "with_metrics/{assay}/{sample}/{quadrant}/{{interval}}.csv"
            ),
                zip,
                assay = ASSAYS,
                sample = SAMPLES,
                quadrant = QUADRANTS
        ),
    output:
        os.path.join(
            config["working_dir"],
            "merged/{interval}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/merge_csvs/{interval}.log"
        ),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 5000
    container:
        config["rocker_tidyverse"]
    script:
        "../scripts/merge_csvs.R"

# Randomise order of videos for 0.5 split into train and test datasets
rule hmm_cocordance_input:
    input:
        rules.merge_csvs.output,
    output:
        os.path.join(
            config["working_dir"],
            "hmm_concordance_in/{interval}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/hmm_cocordance_input/{interval}.log"
        ),
    resources:
        mem_mb = 10000
    container:
        config["rocker_tidyverse"]
    script:
        "../scripts/hmm_concordance_input.R"

rule run_hmm:
    input:
        rules.merge_csvs.output,
    output:
        os.path.join(
            config["working_dir"],
            "hmm_out/{interval}/{variables}/{n_states}.csv"
        ),        
    log:
        os.path.join(
            config["working_dir"],
            "logs/hmm/{interval}/{variables}/{n_states}.log"
        ),
    params:
        n_states = "{n_states}",
        variables = lambda wildcards: config["hmm_variables"][wildcards.variables]
    resources:
        # start at 10000
        mem_mb = lambda wildcards, attempt: attempt * 18000,
    container:
        config["hmmlearn"]
    script:
        "../scripts/run_hmm.py"

rule hmm_concordance:
    input:
        rules.merge_csvs.output,
    output:
        os.path.join(
            config["working_dir"],
            "hmm_concordance/{interval}/{variables}/{n_states}.csv"
        ),        
    log:
        os.path.join(
            config["working_dir"],
            "logs/hmm_concordance_out/{interval}/{variables}/{n_states}.log"
        ),
    params:
        n_states = "{n_states}",
        variables = lambda wildcards: config["hmm_variables"][wildcards.variables]
    resources:
        # start at 10000
        mem_mb = lambda wildcards, attempt: attempt * 18000,
    container:
        config["hmmlearn"]
    script:
        "../scripts/hmm_concordance.py"    
