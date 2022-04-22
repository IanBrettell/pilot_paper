# Recode the states for group B so that they're effectively the same state as group A
rule recode_concordance:
    input:
        rules.hmm_concordance_out.output,
    output:
        rds = os.path.join(
            config["working_dir"],
            "hmm_concordance_recode/{interval}/{variables}/{n_states}.rds"
        ),
        confmat_png = "book/figs/hmm_concordance/{interval}/{variables}/{n_states}/conf_mat.png",
        confmat_pdf = "book/figs/hmm_concordance/{interval}/{variables}/{n_states}/conf_mat.pdf",
        polar_png = "book/figs/hmm_concordance/{interval}/{variables}/{n_states}/eye_facet.png",
        polar_pdf = "book/figs/hmm_concordance/{interval}/{variables}/{n_states}/eye_facet.pdf"
    log:
        os.path.join(
            config["working_dir"],
            "logs/recode_concordance/{interval}/{variables}/{n_states}.log"
        ),
    params:
        interval = "{interval}",
        n_states = "{n_states}",
    resources:
        mem_mb = 20000,
    container:
        config["R"]
    script:
        "../scripts/recode_concordance.R"

rule run_kruskal_wallis:
    input:
        expand(rules.run_hmm.output,
                interval = config["seconds_interval"],
                variables = config["hmm_variables"],
                n_states = config["n_states"]               
        ),
    output:
        os.path.join(
            config["working_dir"],
            "kruskal_wallis/out.rds"
        ),   
    log:
        os.path.join(
            config["working_dir"],
            "logs/run_kruskal_wallis/all.log"
        ),
    resources:
        mem_mb = 30000
    container:
        config["R"]
    script:
        "../scripts/run_kruskal_wallis.R"

rule compare_params:
    input:
        conc = expand(rules.recode_concordance.output.rds,
            interval = config["seconds_interval"],
            variables = config["hmm_variables"],
            n_states = config["n_states"]
        ),
        kw = rules.run_kruskal_wallis.output,
    output:
        png = "book/figs/compare_param_conc.png",
        pdf = "book/figs/concordance_v_kw.pdf"
    log:
        os.path.join(
            config["working_dir"],
            "logs/compare_param_concordance/all.log"
        ),
    resources:
        mem_mb = 20000,
    container:
        config["R"]
    script:
        "../scripts/concordance_v_kw.R"

    
