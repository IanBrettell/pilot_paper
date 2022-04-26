rule hmm_final:
    input:
        rules.run_hmm.output
    output:
        polar_all = "book/figs/hmm_final/{interval}/{variables}/{n_states}.png",
        freq_by_line = "book/figs/hmm_final/{interval}/{variables}/{n_states}.pdf"
    log:
        os.path.join(
            config["working_dir"],
            "logs/hmm_polar_final/{interval}/{variables}/{n_states}.log"
        ),
    resources:
        mem_mb = 30000
    container:
        config["R"]
    script:
        "../scripts/hmm_final.R"