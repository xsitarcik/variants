rule ivar__create_consensus:
    input:
        bam=infer_bam_for_sample_and_ref,
        bai=infer_bai_for_sample_and_ref,
    output:
        consensus=report(
            "results/consensus/ivar/{reference}/{sample}.fa",
            category="{sample} - {reference}",
            labels={
                "Type": "Consensus",
            },
        ),
    params:
        samtools_params=parse_samtools_params_for_consensus(),
        ivar_params=parse_ivar_params_for_consensus(),
    log:
        "logs/ivar/create_consensus/{reference}/{sample}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.12.7/wrappers/ivar/consensus"
