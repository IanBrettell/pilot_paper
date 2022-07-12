# Create plots for covariate effects with p-values
rule covariate_effects:
    input:
        rules.merge_csvs.output,
    output:
        of = "book/figs/covariate_effects/{interval}/{variables}/{n_states}/covariate_effects_of.png",
        no = "book/figs/covariate_effects/{interval}/{variables}/{n_states}/covariate_effects_no.png",
        no_split_by_assay = "book/figs/covariate_effects/{interval}/{variables}/{n_states}/covariate_effects_no-split-by-assay.png",
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

# Create polar plots
rule hmm_final:
    input:
        rules.run_hmm.output
    output:
        polar_all_dge = "book/figs/paper_final/{interval}/{variables}/{n_states}/polar_all_dge.png",
        polar_all_sge = "book/figs/paper_final/{interval}/{variables}/{n_states}/polar_all_sge.png",
        polar_box_dge = "book/figs/paper_final/{interval}/{variables}/{n_states}/polar_box_dge.png",
        polar_box_sge = "book/figs/paper_final/{interval}/{variables}/{n_states}/polar_box_sge.png",
        polar_box_dge_sge = "book/figs/paper_final/{interval}/{variables}/{n_states}/polar_box_dge_sge.png",
        # NOTE: the next four don't work within the script, only in RStudio
        polar_all_dge_sig_of = "book/figs/paper_final/{interval}/{variables}/{n_states}/polar_all_dge_sig_of.png",
        polar_all_dge_sig_no = "book/figs/paper_final/{interval}/{variables}/{n_states}/polar_all_dge_sig_no.png",
        polar_all_sge_sig_of = "book/figs/paper_final/{interval}/{variables}/{n_states}/polar_all_sge_sig_of.png",
        polar_all_sge_sig_no = "book/figs/paper_final/{interval}/{variables}/{n_states}/polar_all_sge_sig_no.png",
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
        mem_mb = 10000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/spatial_dependence.R"
    
rule time_dependence:
    input:
        rules.run_hmm.output,
    output:
        dge = "book/figs/time_dependence/{variables}/{interval}_{n_states}_dge.png",
        sge = "book/figs/time_dependence/{variables}/{interval}_{n_states}_sge.png",
    log:
        os.path.join(
            config["working_dir"],
            "logs/time_dependence/{interval}/{variables}/{n_states}.log"
        ),
    params:
        n_states = "{n_states}",
    resources:
        mem_mb = 10000
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

rule sge_co_occupancy:
    input:
        rules.run_hmm.output,
    output:
        boxplot_all = "book/figs/sge/co-occupancy/{variables}/{interval}_{n_states}_cooc_box_all.png",
        box_and_heat_per_state = "book/figs/sge/co-occupancy/{variables}/{interval}_{n_states}_cooc_box_heat_per-state.png",
    log:
        os.path.join(
            config["working_dir"],
            "logs/sge_co_occupancy/{interval}/{variables}/{n_states}.log"
        ),
    params:
        n_states = "{n_states}",
    resources:
        mem_mb = 3000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/sge_co_occupancy.R"

# Create path videos for both fishes
rule path_frames_to_vid:
    input:
        hmm = rules.run_hmm.output,
        dims = rules.get_split_video_dims.output,
    output:
        os.path.join(
            config["working_dir"],
            "path_vids/{interval}/{variables}/{n_states}/{assay}/{sample}.avi"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/path_frames_to_vid/{interval}/{variables}/{n_states}/{assay}/{sample}.log"
        ),
    params:
        assay = "{assay}",
        sample = "{sample}",
        fps = get_fps,
        tmpdir = config["tmpdir"]
    resources:
        mem_mb = 5000,
    container:
        config["opencv"]
    script:
        "../scripts/path_frames_to_vid.py"

# Creates point/path plots for each ref/test fish
# With each point coloured by HMM state
rule hmm_path_frames_to_vid:
    input:
        hmm = rules.run_hmm.output,
        dims = rules.get_split_video_dims.output,
    output:
        os.path.join(
            config["working_dir"],
            "hmm_path_vids/{interval}/{variables}/{n_states}/{assay}/{sample}_{ref_test}.avi"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/hmm_path_frames_to_vid/{interval}/{variables}/{n_states}/{assay}/{sample}_{ref_test}.log"
        ),
    params:
        assay = "{assay}",
        sample = "{sample}",
        ref_test = "{ref_test}",
        fps = get_fps,
        tmpdir = os.path.join(
            config["working_dir"],
            "tmp_frames/{interval}/{variables}/{n_states}"
        ),
    resources:
        mem_mb = 5000,
    container:
        config["opencv"]
    script:
        "../scripts/hmm_path_frames_to_vid.py"

# Create four-panel videos with labelled videos (TL), path videos (TR),
# and HMM path videos for test (BL) and ref (BR) fishes.
rule compile_four_panel_vid:
    input:
        labels = rules.stitch_tracked_vids.output[0],
        paths = rules.path_frames_to_vid.output[0],
        hmm_ref = os.path.join(
            config["working_dir"],
            "hmm_path_vids/{interval}/{variables}/{n_states}/{assay}/{sample}_ref.avi"
        ),
        hmm_test = os.path.join(
            config["working_dir"],
            "hmm_path_vids/{interval}/{variables}/{n_states}/{assay}/{sample}_test.avi"
        ),
    output:
        os.path.join(
            config["working_dir"],
            "four_panel_vids/{interval}/{variables}/{n_states}/{assay}/{sample}.avi"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/compile_four_panel_vid/{interval}/{variables}/{n_states}/{assay}/{sample}.log"
        ),
    params:
        fps = get_fps,
    resources:
        mem_mb = 5000,
    container:
        config["opencv"]
    script:
        "../scripts/compile_four_panel_vid.py"

rule four_panel_short:
    input:
        rules.compile_four_panel_vid.output,
    output:
        avi = os.path.join(
            config["working_dir"],
            "four_panel_short_vids/{interval}/{variables}/{n_states}/{assay}/{sample}.avi"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/four_panel_short/{interval}/{variables}/{n_states}/{assay}/{sample}.log"
        ),
    params:
        tot_sec = 60,
    resources:
        mem_mb = 5000,
    container:
        config["opencv"]
    script:
        "../scripts/four_panel_short.py"         

# Pull single frame with labeled video on left and paths on right
rule get_labels_and_paths_frame_grab:
    input:
        labels = rules.stitch_tracked_vids.output[0],
        paths = rules.path_frames_to_vid.output[0],
    output:
        fig = "book/figs/labs_and_paths_frame_grabs/{interval}/{variables}/{n_states}/{assay}/{sample}_{second}.png"
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_labels_and_paths_frame_grab/{interval}/{variables}/{n_states}/{assay}/{sample}_{second}.log"
        ),
    params:
        target_second = "{second}"
    resources:
        mem_mb = 500,
    container:
        config["opencv"]
    script:
        "../scripts/get_labels_and_paths_frame_grab.py"

rule compose_setup_figure:
    input:
        frame_grab = rules.get_labels_and_paths_frame_grab.output,
        setup_pic = "book/figs/misc/setup_picture.pdf",
    output:
        fig = "book/figs/misc/setup_fig/{interval}/{variables}/{n_states}/{assay}/{sample}/setup_pic_with_frame_{second}.png",
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_labels_and_paths_frame_grab/{interval}/{variables}/{n_states}/{assay}/{sample}_{second}.log"
        ),
    params:
        setup_png = "book/figs/misc/setup_picture.png",
    resources:
        mem_mb = 5000,
    container:
        config["R_4.2.0"]
    script:
        "../scripts/compose_setup_figure.R"

# Send plots to be used for the paper to the Google Drive
rule send_plots_to_google_drive:
    input:
        setup = expand(rules.compose_setup_figure.output.fig,
                interval = 0.08,
                variables = "dist_angle",
                n_states = 15,
                assay = "open_field",
                sample = "20190613_1054_icab_hdr_R",
                second = 110         
        ),
        covariates = expand(rules.covariate_effects.output.no_split_by_assay,
                interval = 0.08,
                variables = "dist_angle",
                n_states = 15,            
        ),
        compare_params = rules.compare_params.output.png,
        polar = expand(rules.hmm_final.output.polar_all_dge,
                interval = 0.08,
                variables = "dist_angle",
                n_states = 15,        
        ),
        time_dge = expand(rules.time_dependence.output.dge,
                interval = 0.08,
                variables = "dist_angle",
                n_states = 15            
        ),
        time_sge = expand(rules.time_dependence.output.sge,
                interval = 0.08,
                variables = "dist_angle",
                n_states = 15            
        ),        
    output:
        touch(
            os.path.join(
                config["working_dir"],
                "logs/send_plots_to_google_drive/all.done"
            )
        )
    log:
        os.path.join(
            config["working_dir"],
            "logs/send_plots_to_google_drive/all.log"
        ),
    params:
        drive_dir = config["google_drive_dir"]
    conda:
        "../envs/rclone.yaml"
    resources:
        mem_mb = 1000
    shell:
        """
        for i in $(echo {input}); do \
            rclone copy $i {params.drive_dir}/ ; \
        done
        """
