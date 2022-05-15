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
        polar_all = "book/figs/paper_final/{interval}/{variables}/{n_states}/polar_all.png",
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

rule path_videos:
    input:
        hmm = rules.run_hmm.output,
        dims = rules.get_split_video_dims.output,
    output:
        os.path.join(
            config["working_dir"],
            "path_videos/{interval}/{variables}/{n_states}/{assay}/{sample}.avi"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/path_videos/{interval}/{variables}/{n_states}/{assay}/{sample}.log"
        ),
    params:
        assay = "{assay}",
        sample = "{sample}",
        interval = "{interval}"
    resources:
        mem_mb = 200000,
        tmpdir = config["tmpdir"]
    container:
        config["R_4.2.0"]
    script:
        "../scripts/path_videos.R"

rule hmm_path_videos:
    input:
        hmm = rules.run_hmm.output,
        dims = rules.get_split_video_dims.output,
    output:
        os.path.join(
            config["working_dir"],
            "hmm_path_videos/{interval}/{variables}/{n_states}/{assay}/{sample}_{ref_test}.avi"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/hmm_path_videos/{interval}/{variables}/{n_states}/{assay}/{sample}_{ref_test}.log"
        ),
    params:
        assay = "{assay}",
        sample = "{sample}",
        ref_test = "{ref_test}",
        interval = "{interval}"
    resources:
        mem_mb = 80000,
        tmpdir = config["tmpdir"]
    container:
        config["R_4.2.0"]
    script:
        "../scripts/hmm_path_videos.R"

#rule combine_tracked_and_path_vids:
#    input:
#        vid = rules.stitch_tracked_vids.output,
#        path = rules.path_videos.output,
#        hmm = expand(os.path.join(
#            config["working_dir"],
#            "hmm_path_videos/{interval}/{variables}/{n_states}/{assay}/{sample}_{ref_test}.avi"
#            ),
#                ref_test = REF_TEST
#        ),
#    output:
#        os.path.join(
#            config["working_dir"],
#            "combined_videos/{interval}/{variables}/{n_states}/{assay}/{sample}.avi"
#        ),        
#    log:
#        os.path.join(
#            config["working_dir"],
#            "logs/combine_tracked_and_path_vids/{interval}/{variables}/{n_states}/{assay}/{sample}.log"
#        ),
#    resources:
#        mem_mb = 20000
#    container:
#        config["opencv"]
#    script:
#        "../scripts/combine_tracked_and_path_vids.py"    
