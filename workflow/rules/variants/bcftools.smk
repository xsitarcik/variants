rule bcftools__mpileup_bcf:
    input:
        alignments=infer_bam_for_sample_and_ref,
        ref=infer_reference_fasta(),
        index=infer_reference_faidx(),
        bai=infer_bai_for_sample_and_ref,
    output:
        pileup=temp("results/variants_bcftools/{reference}/original/{sample}.pileup.bcf"),
    params:
        uncompressed_bcf=True,
        extra=get_bcftools_mpileup_params(),
    threads: min(config["threads"]["bcftools"], config["max_threads"])
    log:
        "logs/variants_bcftools/mpileup_bcf/{reference}/{sample}.log",
    wrapper:
        "v3.3.3/bio/bcftools/mpileup"


rule bcftools__call_variants:
    input:
        pileup="results/variants_bcftools/{reference}/original/{sample}.pileup.bcf",
    output:
        calls=temp("results/variants_bcftools/{reference}/original/{sample}.bcftools.vcf"),
    params:
        caller="--multiallelic-caller",
        extra=get_bcftools_calling_params(),
    threads: min(config["threads"]["bcftools"], config["max_threads"])
    log:
        "logs/variants_bcftools/call_variants/{reference}/{sample}.log",
    wrapper:
        "v3.3.3/bio/bcftools/call"


rule gatk__prepare_vcf:
    input:
        vcf="results/variants_bcftools/{reference}/original/{sample}.bcftools.vcf",
        ref=infer_reference_fasta(),
        ref_dct=infer_reference_dict(),
    output:
        vcf=temp("results/variants_bcftools/{reference}/original/{sample}.vcf"),
    log:
        "logs/variants_bcftools/prepare_vcf/{reference}/{sample}.log",
    params:
        extra="",  #"--select-type-to-include SNP",  # optional filter arguments, see GATK docs
        java_opts="",  # optional
    wrapper:
        "v3.3.3/bio/gatk/selectvariants"


rule bcftools__normalize_vcf:
    input:
        vcf="results/variants_bcftools/{reference}/original/{sample}.vcf",
        ref=infer_reference_fasta(),
        index=infer_reference_faidx(),
    output:
        temp("results/variants_bcftools/{reference}/normalized/{sample}.norm.bcf"),
    params:
        uncompressed_bcf=True,
        extra=get_bcftools_norm_params(),
    log:
        "logs/variants_bcftools/normalize_vcf/{reference}/{sample}.log",
    threads: min(config["threads"]["bcftools"], config["max_threads"])
    wrapper:
        "v3.3.3/bio/bcftools/norm"


rule bcftools__annotate_vcf:
    input:
        "results/variants_bcftools/{reference}/normalized/{sample}.norm.bcf",
    output:
        temp("results/variants_bcftools/{reference}/normalized/{sample}.annot.vcf"),
    log:
        "logs/variants_bcftools/annotate_vcf/{reference}/{sample}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -O v -o {output} {input} 1> {log} 2>&1"


rule bcftools__fill_tags:
    input:
        "results/variants_bcftools/{reference}/normalized/{sample}.annot.vcf",
    output:
        temp("results/variants_bcftools/{reference}/normalized/{sample}.vcf"),
    log:
        "logs/variants_bcftools/tag_fill/{reference}/{sample}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools plugin fill-tags {input} -O v -o {output} -- --tags 'AF,AC,AN' 1> {log} 2>&1"


rule bcftools__filter_vcf:
    input:
        "results/variants_bcftools/{reference}/normalized/{sample}.vcf",
    output:
        temp("results/variants_bcftools/{reference}/filtered/{sample}.filtered.bcf"),
    params:
        extra=get_bcftools_filter_params(),
    log:
        "logs/variants_bcftools/filter_vcf/{reference}/{sample}.log",
    threads: min(config["threads"]["bcftools"], config["max_threads"])
    wrapper:
        "v3.3.3/bio/bcftools/filter"


rule bcftools__view_filtered_vcf:
    input:
        "results/variants_bcftools/filtered/{reference}/{sample}.filtered.bcf",
    output:
        report(
            "results/variants_bcftools/filtered/{reference}/{sample}.vcf",
            category="{sample} - {reference}",
            labels={
                "Type": "Variants - bcftools/filtered",
            },
        ),
    params:
        extra="--trim-alt-alleles" if config["variants__postprocess"]["trim_alt_alleles"] else "",
    log:
        "logs/variants_bcftools/view_filtered_vcf/{reference}/{sample}.log",
    threads: min(config["threads"]["bcftools"], config["max_threads"])
    wrapper:
        "v3.3.3/bio/bcftools/view"
