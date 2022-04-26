rule hmm_final:
    input:
        rules.run_hmm.output
    output:
        polar_box_dge = "book/figs/paper_final/{interval}/{variables}/{n_states}/polar_box_dge.png",
        polar_box_sge = "book/figs/paper_final/{interval}/{variables}/{n_states}/polar_box_sge.png",
        tile_dge = "book/figs/paper_final/{interval}/{variables}/{n_states}/tile_dge.png",
        tile_sge = "book/figs/paper_final/{interval}/{variables}/{n_states}/tile_sge.png",
    log:
        os.path.join(
            config["working_dir"],
            "logs/hmm_final/{interval}/{variables}/{n_states}.log"
        ),
    params:
        interval = "{interval}",
        n_states = "{n_states}",
    resources:
        mem_mb = 30000
    container:
        config["R"]
    script:
        "../scripts/hmm_final.R"