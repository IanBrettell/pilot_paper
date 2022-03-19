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
        # start at 5000
        mem_mb = lambda wildcards, attempt: attempt * 5000,
    container:
        config["hmmlearn"]
    script:
        "../scripts/run_hmm.py"

# We also want to test the concordance between the assigned states
# when the HMM is trained on a different dataset
# Randomise order of videos for 0.5 split into train and test datasets
rule hmm_concordance_in:
    input:
        rules.merge_csvs.output,
    output:
        A = os.path.join(
            config["working_dir"],
            "hmm_concordance_in/{interval}/A.csv"
        ),
        B = os.path.join(
            config["working_dir"],
            "hmm_concordance_in/{interval}/B.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/hmm_cocordance_input/{interval}.log"
        ),
    resources:
        mem_mb = 5000
    container:
        config["rocker_tidyverse"]
    script:
        "../scripts/hmm_concordance_input.R"

# Run concordance
rule hmm_concordance_out:
    input:
        A = rules.hmm_concordance_in.output.A,
        B = rules.hmm_concordance_in.output.B,
    output:
        os.path.join(
            config["working_dir"],
            "hmm_concordance_out/{interval}/{variables}/{n_states}.csv"
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
        # start at 5000
        mem_mb = lambda wildcards, attempt: attempt * 10000,
    container:
        config["hmmlearn"]
    script:
        "../scripts/hmm_concordance.py"    
