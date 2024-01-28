rule ivar__get_variants:
    input:
        bam="results/mapping/deduplicated/{sample}.bam",
        bai="results/mapping/deduplicated/{sample}.bam.bai",
        ref=infer_reference_fasta,
    output:
        tsv=temp("results/variants_ivar/{reference}/{sample}/all.tsv"),
    params:
        samtools_params=parse_samtools_params_for_variants(),
        ivar_params=parse_ivar_params_for_variants(),
    log:
        "logs/ivar/get_variants/{reference}/{sample}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.12.7/wrappers/ivar/variants"


rule ivar__variants_to_vcf:
    input:
        "results/variants_ivar/{reference}/{sample}/all.tsv",
    output:
        all=report(
            "results/variants_ivar/{reference}/{sample}/all.vcf",
            category="{sample} - {reference}",
            labels={
                "Type": "Variants - Ivar vcf",
            },
        ),
        filtered="results/variants_ivar/{reference}/{sample}/passed_only.vcf",
    log:
        "logs/ivar/variants_to_vcf/{reference}/{sample}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.12.7/wrappers/ivar/vcf_convert"


rule custom__compute_mixed_positions:
    input:
        "results/variants_ivar/{reference}/{sample}/all.tsv",
    output:
        mixed_positions="results/variants_ivar/{reference}/{sample}/filtered.tsv",
        readcount=temp("results/variants_ivar/{reference}/{sample}/filtered_count.tsv"),
    params:
        unpack(get_ivar_postfilter_params()),
    log:
        "logs/mixed_positions/{reference}/{sample}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.12.7/wrappers/ivar/mixed_positions"


rule report__ivar_variants_to_html:
    input:
        "results/variants_ivar/{reference}/{sample}/all.tsv",
    output:
        report(
            "results/variants_ivar/{reference}/{sample}/all.html",
            category="{sample} - {reference}",
            labels={
                "Type": "Variants - Ivar table",
            },
        ),
    log:
        "logs/report/ivar_variants/{reference}/{sample}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.12.7/wrappers/ivar/html_convert"


rule report__mixed_positions_to_html:
    input:
        "results/variants_ivar/{reference}/{sample}/filtered.tsv",
    output:
        report(
            "results/variants_ivar/{reference}/{sample}/filtered.html",
            category="{sample} - {reference}",
            labels={
                "Type": "Filtered variants - Ivar",
            },
        ),
    log:
        "logs/report/mixed_positions/{reference}/{sample}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.12.7/wrappers/ivar/html_convert"
