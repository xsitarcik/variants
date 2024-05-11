rule bcftools__stats:
    input:
        vcf="results/variants/{reference}/{sample}/{tool}_{step}.vcf",
        ref=infer_reference_fasta,
    output:
        report(
            "results/variants/{reference}/{sample}/stats/{tool}_{step}.txt",
            category="Variants - {reference} - stats",
            labels={
                "Sample": "{sample}",
                "Type": "Stats for {tool} - {step}",
            },
        ),
    log:
        "logs/variants_stats/{reference}/{tool}/{sample}_{step}.log",
    wrapper:
        "v3.10.2/bio/bcftools/stats"
