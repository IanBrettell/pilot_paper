rule sge_co_occupancy:
    input:
        rules.run_hmm.output,
    output:
        heatmaps = "book/figs/sge/co-occupancy/{interval}/{variables}/{n_states}/cooc_heat.png",
        boxplot_all = "book/figs/sge/co-occupancy/{interval}/{variables}/{n_states}/cooc_box_all.png",
        boxplot_per_state = "book/figs/sge/co-occupancy/{interval}/{variables}/{n_states}/cooc_box_per-state.png",
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