rule ivar__get_variants:
    input:
        bam="results/mapping/deduplicated/{sample}.bam",
        bai="results/mapping/deduplicated/{sample}.bam.bai",
        ref=get_reference_fasta(),
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
            category="{sample}",
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
        "results/variants_ivar/{sample}/all.tsv",
    output:
        mixed_positions="results/variants_ivar/{sample}/mixed_positions.tsv",
        readcount=temp("results/variants_ivar/{sample}/mixed_positions_count.tsv"),
    params:
        alt_depth=config["mixed_positions_params"]["filtering"]["min_alt_depth"],
        min_alt_freq=config["mixed_positions_params"]["filtering"]["min_alt_freq"],
        max_alt_freq=config["mixed_positions_params"]["filtering"]["max_alt_freq"],
        total_depth=config["mixed_positions_params"]["filtering"]["min_total_depth"],
    log:
        "logs/mixed_positions/{sample}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.12.7/wrappers/ivar/mixed_positions"


rule report__ivar_variants_to_html:
    input:
        "results/variants_ivar/{sample}/all.tsv",
    output:
        report(
            "results/variants_ivar/{sample}/all.html",
            category="{sample}",
            labels={
                "Type": "Variants - Ivar table",
            },
        ),
    log:
        "logs/report/ivar_variants/{sample}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.12.7/wrappers/ivar/html_convert"


rule report__mixed_positions_to_html:
    input:
        "results/variants_ivar/{sample}/mixed_positions.tsv",
    output:
        report(
            "results/variants_ivar/{sample}/mixed_positions.html",
            category="{sample}",
            labels={
                "Type": "Mixed positions - Ivar",
            },
        ),
    log:
        "logs/report/mixed_positions/{sample}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.12.7/wrappers/ivar/html_convert"
