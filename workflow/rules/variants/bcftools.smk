rule bcftools__mpileup_bcf:
    input:
        alignments=infer_bam_for_sample_and_ref,
        ref=infer_reference_fasta,
        index=infer_reference_faidx,
        bai=infer_bai_for_sample_and_ref,
    output:
        pileup=temp("results/variants/{reference}/{sample}/bcftools_all.pileup.bcf"),
    params:
        uncompressed_bcf=True,
        extra=get_bcftools_mpileup_params(),
    threads: get_threads_for_bcftools()
    log:
        "logs/variants_bcftools/mpileup_bcf/{reference}/{sample}.log",
    wrapper:
        "v3.10.2/bio/bcftools/mpileup"


rule bcftools__call_variants:
    input:
        pileup="results/variants/{reference}/{sample}/bcftools_all.pileup.bcf",
    output:
        calls=temp("results/variants/{reference}/{sample}/bcftools_all.bcftools.vcf"),
    params:
        caller="--multiallelic-caller",
        extra=get_bcftools_calling_params(),
    threads: get_threads_for_bcftools()
    log:
        "logs/variants_bcftools/call_variants/{reference}/{sample}.log",
    wrapper:
        "v3.10.2/bio/bcftools/call"


rule bcftools__gatk_prepare_vcf:
    input:
        vcf="results/variants/{reference}/{sample}/bcftools_all.bcftools.vcf",
        ref=infer_reference_fasta,
        ref_dct=infer_reference_dict,
    output:
        vcf=temp("results/variants/{reference}/{sample}/bcftools_all.vcf"),
        idx=temp("results/variants/{reference}/{sample}/bcftools_all.vcf.idx"),
    log:
        "logs/variants_bcftools/gatk_prepare_vcf/{reference}/{sample}.log",
    params:
        extra="",
        java_opts="",
    wrapper:
        "v3.10.2/bio/gatk/selectvariants"


rule bcftools__normalize_vcf:
    input:
        vcf="results/variants/{reference}/{sample}/bcftools_all.vcf",
        ref=infer_reference_fasta,
        index=infer_reference_faidx,
    output:
        temp("results/variants/{reference}/{sample}/bcftools_norm.bcf"),
    params:
        uncompressed_bcf=True,
        extra=get_bcftools_norm_params("bcftools"),
    log:
        "logs/variants_bcftools/normalize_vcf/{reference}/{sample}.log",
    threads: get_threads_for_bcftools()
    wrapper:
        "v3.10.2/bio/bcftools/norm"


rule bcftools__annotate_vcf:
    input:
        "results/variants/{reference}/{sample}/bcftools_norm.bcf",
    output:
        temp("results/variants/{reference}/{sample}/bcftools_annot.vcf"),
    log:
        "logs/variants_bcftools/annotate_vcf/{reference}/{sample}.log",
    conda:
        "../../envs/bcftools.yaml"
    shell:
        "bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -O v -o {output} {input} 1> {log} 2>&1"


rule bcftools__fill_tags:
    input:
        "results/variants/{reference}/{sample}/bcftools_annot.vcf",
    output:
        temp("results/variants/{reference}/{sample}/bcftools_filltags.vcf"),
    log:
        "logs/variants_bcftools/fill_tags/{reference}/{sample}.log",
    params:
        tags=parse_tags_for_bcftools_fill_tags("bcftools"),
    conda:
        "../../envs/bcftools.yaml"
    shell:
        "bcftools plugin fill-tags {input} -O v -o {output} -- {params.tags} 1> {log} 2>&1"


rule bcftools__filter_vcf:
    input:
        "results/variants/{reference}/{sample}/bcftools_filltags.vcf",
    output:
        temp("results/variants/{reference}/{sample}/bcftools_filtered.bcf"),
    params:
        extra=get_bcftools_filter_params("bcftools"),
    log:
        "logs/variants_bcftools/filter_vcf/{reference}/{sample}.log",
    threads: get_threads_for_bcftools()
    wrapper:
        "v3.10.2/bio/bcftools/filter"


rule bcftools__view_filtered:
    input:
        "results/variants/{reference}/{sample}/bcftools_filtered.bcf",
    output:
        report(
            "results/variants/{reference}/{sample}/bcftools_filtered.vcf",
            category="Variants - {reference}",
            labels={
                "Sample": "{sample}",
                "Type": "BCFtools - filtered",
            },
        ),
    params:
        extra=get_bcftools_view_filter_params("bcftools"),
    log:
        "logs/variants_bcftools/view_filtered_vcf/{reference}/{sample}.log",
    threads: get_threads_for_bcftools()
    wrapper:
        "v3.10.2/bio/bcftools/view"
