# https://snakemake-wrappers.readthedocs.io/en/stable/meta-wrappers/gatk_mutect2_calling.html


rule picard__replace_read_groups:
    input:
        bam=infer_bam_for_sample_and_ref,
    output:
        temp("results/variants_mutect2/{reference}/{sample}/fixed.bam"),
    threads: 1
    resources:
        mem_mb=1024,
    log:
        "logs/picard/replace_rg/{sample}.log",
    params:
        # Required for GATK
        extra="--RGLB lib1 --RGPL illumina --RGPU {sample} --RGSM {sample}",
    wrapper:
        "v3.8.0/bio/picard/addorreplacereadgroups"


rule sambamba__index_bam:
    input:
        "results/variants_mutect2/{reference}/{sample}/fixed.bam",
    output:
        "results/variants_mutect2/{reference}/{sample}/fixed.bam.bai",
    threads: 1
    log:
        "logs/sambamba/index/{sample}.log",
    params:
        extra="",
    wrapper:
        "v3.8.0/bio/sambamba/index"


rule mutect2_call:
    input:
        fasta=infer_reference_fasta,
        fasta_dict=infer_reference_dict,
        fasta_fai=infer_reference_faidx,
        map="results/variants_mutect2/{reference}/{sample}/fixed.bam",
        map_idx="results/variants_mutect2/{reference}/{sample}/fixed.bam.bai",
    output:
        vcf="results/variants_mutect2/{reference}/{sample}/original.vcf",
        tbi="results/variants_mutect2/{reference}/{sample}/original.vcf.tbi",
        f1r2="results/variants_mutect2/{reference}/{sample}/counts.f1r2.tar.gz",
    threads: 1
    resources:
        mem_mb=1024,
    params:
        extra="--create-output-variant-index --tumor-sample {sample} ",  # TODO
    log:
        "logs/mutect/{sample}.log",
    wrapper:
        "v3.8.0/bio/gatk/mutect"


rule gatk_get_pileup_summaries:
    input:
        bam="results/variants_mutect2/{reference}/{sample}/fixed.bam",
        bai_bai="results/variants_mutect2/{reference}/{sample}/fixed.bam.bai",
        variants="results/variants_mutect2/{reference}/{sample}/original.vcf",
        variants_tbi="results/variants_mutect2/{reference}/{sample}/original.vcf.tbi",
    output:
        "results/variants_mutect2/{reference}/{sample}/summary.table",
    threads: 1
    resources:
        mem_mb=1024,
    params:
        extra="",
    log:
        "logs/summary/{sample}.log",
    wrapper:
        "v3.8.0/bio/gatk/getpileupsummaries"


rule gatk_calculate_contamination:
    input:
        tumor="results/variants_mutect2/{reference}/{sample}/summary.table",
    output:
        temp("results/variants_mutect2/{reference}/{sample}/summary.pileups.table"),
    threads: 1
    resources:
        mem_mb=1024,
    log:
        "logs/contamination/{sample}.log",
    params:
        extra="",
    wrapper:
        "v3.8.0/bio/gatk/calculatecontamination"


rule gatk_learn_read_orientation_model:
    input:
        f1r2="results/variants_mutect2/{reference}/{sample}/counts.f1r2.tar.gz",
    output:
        temp("results/variants_mutect2/{reference}/{sample}/artifacts_prior.tar.gz"),
    threads: 1
    resources:
        mem_mb=1024,
    params:
        extra="",
    log:
        "logs/learnreadorientationbias/{sample}.log",
    wrapper:
        "v3.8.0/bio/gatk/learnreadorientationmodel"


rule filter_mutect_calls:
    input:
        vcf="results/variants_mutect2/{reference}/{sample}/original.vcf",
        ref=infer_reference_fasta,
        ref_dict=infer_reference_dict,
        ref_fai=infer_reference_faidx,
        bam="results/variants_mutect2/{reference}/{sample}/fixed.bam",
        bam_bai="results/variants_mutect2/{reference}/{sample}/fixed.bam.bai",
        contamination="results/variants_mutect2/{reference}/{sample}/summary.pileups.table",
        f1r2="results/variants_mutect2/{reference}/{sample}/artifacts_prior.tar.gz",
    output:
        vcf="results/variants_mutect2/{reference}/{sample}/filtered.vcf",
        vcf_idx="results/variants_mutect2/{reference}/{sample}/filtered.vcf.tbi",
    threads: 1
    resources:
        mem_mb=1024,
    log:
        "logs/gatk/filter/{sample}.log",
    params:
        extra="--microbial-mode --create-output-variant-index --min-median-mapping-quality 35 --max-alt-allele-count 3",
        java_opts="",
    wrapper:
        "v3.8.0/bio/gatk/filtermutectcalls"
