rule covariate_effects:
    input:
        rules.merge_csvs.output,
    output:
        fig = "book/figs/covariate_effects/{interval}/{variables}/{n_states}/covariate_effects.png",
    log:
        os.path.join(
            config["working_dir"],
            "logs/covariate_effects/{interval}/{variables}/{n_states}.log"
        ),
    resources:
        mem_mb = 5000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/covariate_effects.R"

rule hmm_final:
    input:
        rules.run_hmm.output
    output:
        polar_box_dge = "book/figs/paper_final/{interval}/{variables}/{n_states}/polar_box_dge.png",
        polar_box_sge = "book/figs/paper_final/{interval}/{variables}/{n_states}/polar_box_sge.png",
        polar_box_dge_sge = "book/figs/paper_final/{interval}/{variables}/{n_states}/polar_box_dge_sge.png",
        tile_dge = "book/figs/paper_final/{interval}/{variables}/{n_states}/tile_dge.png",
        tile_sge = "book/figs/paper_final/{interval}/{variables}/{n_states}/tile_sge.png",
        tile_dge_sge = "book/figs/paper_final/{interval}/{variables}/{n_states}/tile_dge_sge.png",
    log:
        os.path.join(
            config["working_dir"],
            "logs/hmm_final/{interval}/{variables}/{n_states}.log"
        ),
    params:
        interval = "{interval}",
        n_states = "{n_states}",
    resources:
        mem_mb = 10000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/hmm_final.R"

rule spatial_dependence:
    input:
        rules.run_hmm.output,
    output:
        fig = "book/figs/spatial_dependence/{variables}/{interval}_{n_states}_spatial_dependence.png"
    log:
        os.path.join(
            config["working_dir"],
            "logs/spatial_dependence/{interval}/{variables}/{n_states}.log"
        ),
    params:
        n_states = "{n_states}",
    resources:
        mem_mb = 3000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/spatial_dependence.R"
    
rule time_dependence:
    input:
        rules.run_hmm.output,
    output:
        fig = "book/figs/time_dependence/{variables}/{interval}_{n_states}_time_dependence.png",
    log:
        os.path.join(
            config["working_dir"],
            "logs/time_dependence/{interval}/{variables}/{n_states}.log"
        ),
    params:
        n_states = "{n_states}",
    resources:
        mem_mb = 5000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/time_dependence.R"

rule path_plots:
    input:
        data = rules.run_hmm.output,
        dims = rules.get_split_video_dims.output,
        screenshots = expand("results/split_coord_images/{assay}/{{sample}}.png",
            assay = list(set(ASSAYS))
        ),
    output:
        paths = "book/figs/path_plots/{interval}/{variables}/{n_states}/{sample}.png",
    log:
        os.path.join(
            config["working_dir"],
            "logs/path_plots/{interval}/{variables}/{n_states}/{sample}.log"
        ),
    params:
        sample = "{sample}"
    resources:
        mem_mb = 3000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/path_plots.R"

rule tracking_success_plot:
    input:
        rules.tracking_success.output,
    output:
        obj = "results/figs/tracking_success/tracking_success_gg.rds",
        fig = "book/figs/tracking_success/tracking_success.png",
    log:
        os.path.join(
            config["working_dir"],
            "logs/tracking_success_plot/all.log"
        ),
    params:
        sample = "{sample}"
    resources:
        mem_mb = 3000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/tracking_success_plot.R"