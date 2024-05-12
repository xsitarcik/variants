rule bedtools__mask_bam_low_coverage:
    input:
        bam=infer_bam_for_sample_and_ref,
        bai=infer_bai_for_sample_and_ref,
    output:
        bed=temp("results/variants/{reference}/{sample}/low_coverage.bed"),
    log:
        "logs/consensus_bcftools/mask_low_coverage/{reference}/{sample}.log",
    params:
        min_coverage=get_bcftools_consensus_mask(),
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools genomecov -ibam {input.bam} -bga | awk '$4<{params.min_coverage}' > {output.bed} 2>{log}"


rule bcftools__compress_index:
    input:
        vcf=infer_final_vcf,
    output:
        vcf_gz=temp("results/variants/{reference}/{sample}/{tool}_for_consensus.vcf.gz"),
        csi=temp("results/variants/{reference}/{sample}/{tool}_for_consensus.vcf.gz.csi"),
    log:
        "logs/consensus_bcftools/compress_index/{reference}/{sample}_{tool}.log",
    conda:
        "../../envs/bcftools.yaml"
    shell:
        "bcftools view {input} -O z -o {output.vcf_gz} -W > {log} 2>&1"


rule bcftools__consensus:
    input:
        ref=infer_reference_fasta,
        index=infer_reference_faidx,
        vcf_gz="results/variants/{reference}/{sample}/{tool}_for_consensus.vcf.gz",
        csi="results/variants/{reference}/{sample}/{tool}_for_consensus.vcf.gz.csi",
        bed="results/variants/{reference}/{sample}/low_coverage.bed",
    output:
        report(
            "results/consensus/{reference}/{sample}/{tool}.fa",
            category="Consensus - {reference}",
            labels={
                "Sample": "{sample}",
                "Type": "{tool}",
            },
        ),
    log:
        "logs/consensus_bcftools/{reference}/{sample}_{tool}.log",
    wildcard_constraints:
        tool="(freebayes|mutect2|bcftools)",
    params:
        extra=get_bcftools_consensus_extra(),
        prefix=lambda wildcards: f"{wildcards.sample}_",
    conda:
        "../../envs/bcftools.yaml"
    shell:
        "bcftools consensus --prefix {params.prefix} --fasta-ref {input.ref} --mask {input.bed} {params.extra} {input.vcf_gz} 1> {output} 2>{log}"
