rule freebayes__call_vcf:
    input:
        alns=infer_bam_for_sample_and_ref,
        idxs=infer_bai_for_sample_and_ref,
        ref=infer_reference_fasta,
        index=infer_reference_faidx,
    output:
        vcf=temp("results/variants_freebayes/{reference}/{sample}/original.vcf"),
    log:
        "logs/freebayes/call_vcf/{reference}/{sample}.log",
    params:
        normalize=get_bcftools_norm_params("freebayes"),
        extra=parse_freebayes_params(),
    threads: get_threads_for_freebayes()
    resources:
        mem_mb=get_mem_mb_for_freebayes,
    wrapper:
        "v3.7.0/bio/freebayes"  # THIS should stay frozen as v3.8.0 version has memory issues


rule freebayes__annotate_vcf:
    input:
        "results/variants_freebayes/{reference}/{sample}/original.vcf",
    output:
        temp("results/variants_freebayes/{reference}/{sample}/annot.vcf"),
    log:
        "logs/variants_freebayes/annotate_vcf/{reference}/{sample}.log",
    conda:
        "../../envs/bcftools.yaml"
    shell:
        "bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -O v -o {output} {input} 1> {log} 2>&1"


rule freebayes__fill_tags:
    input:
        "results/variants_freebayes/{reference}/{sample}/annot.vcf",
    output:
        temp("results/variants_freebayes/{reference}/{sample}/filled.vcf"),
    log:
        "logs/variants_freebayes/tag_fill/{reference}/{sample}.log",
    conda:
        "../../envs/bcftools.yaml"
    shell:
        "bcftools plugin fill-tags {input} -O v -o {output} -- --tags 'AF,AC,AN' 1> {log} 2>&1"


rule freebayes__filter_vcf:
    input:
        "results/variants_freebayes/{reference}/{sample}/filled.vcf",
    output:
        temp("results/variants_freebayes/{reference}/{sample}/filtered.bcf"),
    params:
        extra=get_bcftools_filter_params("freebayes"),
    log:
        "logs/variants_freebayes/filter_vcf/{reference}/{sample}.log",
    threads: get_threads_for_freebayes()
    wrapper:
        "v3.7.0/bio/bcftools/filter"


rule freebayes__view_filtered_vcf:
    input:
        "results/variants_freebayes/{reference}/{sample}/filtered.bcf",
    output:
        report(
            "results/variants_freebayes/{reference}/{sample}/filtered.vcf",
            category="{sample} - {reference}",
            labels={
                "Type": "Variants - freebayes/filtered",
            },
        ),
    params:
        extra=get_bcftools_view_filter_params("freebayes"),
    log:
        "logs/variants_freebayes/view_filtered_vcf/{reference}/{sample}.log",
    threads: get_threads_for_freebayes()
    wrapper:
        "v3.7.0/bio/bcftools/view"
