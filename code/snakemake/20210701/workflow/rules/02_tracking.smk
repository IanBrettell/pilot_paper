#rule copy_results_from_gcs:
#    output:
#        os.path.join(config["lts_dir"], "videos/tracking_results/_session_{sample}_{quadrant}_{assay}")
#    params:
#        out_dir = lambda wildcards, output: os.path.dirname(output[0])
#    conda:
#        "../envs/gsutil.yaml"
#    threads:
#        4
#    shell:
#        """
#        gsutil -m cp -r {config[gcs_dir]}/results/* {params.out_dir}
#        """