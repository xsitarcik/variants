# https://snakemake-wrappers.readthedocs.io/en/stable/meta-wrappers/gatk_mutect2_calling.html


rule custom__create_wgs_bed_file:
    input:
        fai=infer_reference_faidx,
    output:
        bed="results/variants_mutect2/{reference}/regions.bed",
    log:
        "logs/custom/create_wgs_bed_file/{reference}.log",
    localrule: True
    conda:
        "../../envs/awk_sed.yaml"
    shell:
        "(awk '{{print $1, 1, $2, $1}}' {input.fai} | sed 's/ /\t/g' > {output.bed}) 2> {log}"


rule picard__bed_to_interval_list:
    input:
        bed="results/variants_mutect2/{reference}/regions.bed",
        dict=infer_reference_dict,
    output:
        intervals="results/variants_mutect2/{reference}/regions.interval_list",
    log:
        "logs/picard/bed_to_interval_list/{reference}.log",
    params:
        extra="--SORT true --UNIQUE true",
    resources:
        mem_mb=1024,
    wrapper:
        "v3.8.0/bio/picard/bedtointervallist"


rule picard__replace_read_groups:
    input:
        bam=infer_bam_for_sample_and_ref,
    output:
        temp("results/variants_mutect2/{reference}/{sample}/fixed.bam"),
    threads: get_threads_for_mutect2()
    resources:
        mem_mb=1024,
    log:
        "logs/picard/replace_rg/{reference}/{sample}.log",
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
    threads: get_threads_for_mutect2()
    log:
        "logs/sambamba/index/{reference}/{sample}.log",
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
        intervals="results/variants_mutect2/{reference}/regions.interval_list",
    output:
        vcf="results/variants_mutect2/{reference}/{sample}/original.vcf",
        tbi="results/variants_mutect2/{reference}/{sample}/original.vcf.idx",
        f1r2="results/variants_mutect2/{reference}/{sample}/counts.f1r2.tar.gz",
    threads: get_threads_for_mutect2()
    resources:
        mem_mb=1024,
    params:
        extra="--create-output-variant-index --tumor-sample {sample} ",  # TODO
    log:
        "logs/mutect/{reference}/{sample}.log",
    wrapper:
        "v3.8.0/bio/gatk/mutect"


rule gatk_learn_read_orientation_model:
    input:
        f1r2="results/variants_mutect2/{reference}/{sample}/counts.f1r2.tar.gz",
    output:
        temp("results/variants_mutect2/{reference}/{sample}/artifacts_prior.tar.gz"),
    threads: get_threads_for_mutect2()
    resources:
        mem_mb=1024,
    params:
        extra="",
    log:
        "logs/learnreadorientationbias/{reference}/{sample}.log",
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
        f1r2="results/variants_mutect2/{reference}/{sample}/artifacts_prior.tar.gz",
    output:
        vcf="results/variants_mutect2/{reference}/{sample}/filtered.vcf",
        stats="results/variants_mutect2/{reference}/{sample}/filtered.stats",
    threads: get_threads_for_mutect2()
    resources:
        mem_mb=1024,
    log:
        "logs/gatk/filter/{reference}/{sample}.log",
    params:
        extra="--microbial-mode --create-output-variant-index --min-median-mapping-quality 35 --max-alt-allele-count 3",
        java_opts="",
    wrapper:
        "v3.8.0/bio/gatk/filtermutectcalls"
