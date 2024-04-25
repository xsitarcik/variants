rule multiqc__report:
    input:
        **get_multiqc_inputs(),
    output:
        "results/_aggregation/multiqc.html",
    params:
        use_input_files_only=True,
        extra=f"--config {workflow.basedir}/resources/multiqc.yaml",
    log:
        "logs/multiqc/all.log",
    wrapper:
        "v3.8.0/bio/multiqc"
