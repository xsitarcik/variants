rule multiqc__report:
    input:
        **get_multiqc_inputs(),
        config=f"{workflow.basedir}/resources/multiqc.yaml",
    output:
        "results/_aggregation/multiqc.html",
    params:
        use_input_files_only=True,
        extra=f"-d -dd 1",
    log:
        "logs/multiqc/all.log",
    wrapper:
        "v3.8.0/bio/multiqc"
