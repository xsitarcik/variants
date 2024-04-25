rule multiqc__report:
    input:
        **get_multiqc_inputs_no_variants(),
        config=f"{workflow.basedir}/resources/multiqc.yaml",
    output:
        "results/_aggregation/multiqc.html",
    params:
        use_input_files_only=True,
        extra="",
    log:
        "logs/multiqc/all.log",
    wrapper:
        "v3.8.0/bio/multiqc"


rule multiqc__report_variants:
    input:
        **get_multiqc_inputs_variants(),
        config=f"{workflow.basedir}/resources/multiqc.yaml",
    output:
        "results/_aggregation/multiqc_variants.html",
    params:
        use_input_files_only=True,
        extra=f"-d -dd 2",
    log:
        "logs/multiqc/variants.log",
    wrapper:
        "v3.8.0/bio/multiqc"
