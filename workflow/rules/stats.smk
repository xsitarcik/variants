rule bcftools__stats:
    input:
        vcfs=get_final_vcf_files(),
        ref=infer_reference_fasta,
    output:
        report(
            "results/variants/{reference}/{sample}/stats.txt",
            category="{sample} - {reference}",
            labels={
                "Type": "Variants statistics for final step",
            },
        ),
    log:
        "logs/variants_stats/{reference}/{sample}.log",
    params:
        extra="",
    conda:
        "../../envs/bcftools.yaml"
    shell:
        "bcftools stats --fasta-ref {input.ref} {input.vcfs} > {output} 2> {log}"
