rule freebayes__call_vcf:
    input:
        alns=infer_bam_for_sample_and_ref,
        idxs=infer_bai_for_sample_and_ref,
        ref=infer_reference_fasta,
        index=infer_reference_faidx,
    output:
        vcf="results/variants_freebayes/{reference}/{sample}.vcf",
    log:
        "logs/freebayes/call_vcf/{reference}/{sample}.log",
    params:
        normalize="-a",  # (one of `-a`, `-f`, `-m`, `-D` or `-d` must be used)
        extra=parse_freebayes_params(),
    threads: get_threads_for_freebayes()
    resources:
        mem_mb=get_mem_mb_for_freebayes,
    wrapper:
        "v3.8.0/bio/freebayes"
