rule bcftools__stats:
    input:
        vcf="results/variants/{reference}/{tool}/{sample}_{step}.vcf",
        ref=infer_reference_fasta,
    output:
        report(
            "results/variants/{reference}/{tool}/stats/{sample}_{step}.txt",
            category="{sample} - {reference}",
            labels={
                "Type": "{tool} - {step}",
            },
        ),
    log:
        "logs/variants_stats/{reference}/{tool}/{sample}_{step}.log",
    params:
        extra="",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools stats --fasta-ref {input.ref} {input.vcf} > {output} 2> {log}"
