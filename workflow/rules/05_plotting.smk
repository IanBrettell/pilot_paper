# Recode the states for group B so that they're effectively the same state as group A
rule recode_concordance:
    input:
        rules.hmm_concordance_out.output,
    output:
        csv = os.path.join(
            config["working_dir"],
            "hmm_concordance_recode/{interval}/{variables}/{n_states}.csv"
        ),
        png = "book/figs/hmm_concordance/{interval}/{variables}/{n_states}/eye_facet.png",
        pdf = "book/figs/hmm_concordance/{interval}/{variables}/{n_states}/eye_facet.pdf"
    log:
        os.path.join(
            config["working_dir"],
            "logs/recode_concordance/{interval}/{variables}/{n_states}.log"
        ),
    params:
        n_states = "{n_states}"
    resources:
        mem_mb = 20000,
    container:
        config["R"]
    script:
        "../scripts/recode_concordance.R"

rule concordance_v_kw:
    input:
        expand(rules.recode_concordance.output.csv,
            interval = config["seconds_interval"],
            variables = config["hmm_variables"],
            n_states = config["n_states"]
        ),
    output:
        png = "book/figs/concordance_v_kw.png",
        pdf = "book/figs/concordance_v_kw.pdf"
    log:
        os.path.join(
            config["working_dir"],
            "logs/concordance_v_kw/all.log"
        ),
    resources:
        mem_mb = 20000,
    container:
        config["R"]
    script:
        "../scripts/concordance_v_kw.R"

    
