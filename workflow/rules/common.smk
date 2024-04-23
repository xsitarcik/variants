from snakemake.utils import validate


configfile: "config/config.yaml"


validate(config, "../schemas/config.schema.yaml", set_default=False)


### Layer for adapting other workflows  ###############################################################################


def get_sample_names():
    return mapping_workflow.get_sample_names()


def get_reference_names():
    return mapping_workflow.get_reference_names()


def infer_bam_for_sample_and_ref(wildcards):
    return mapping_workflow.get_input_bam_for_sample_and_ref(wildcards.sample, wildcards.reference)


def infer_bai_for_sample_and_ref(wildcards):
    return mapping_workflow.get_input_bai_for_sample_and_ref(wildcards.sample, wildcards.reference)


def get_reference_dir_for_name(reference: str):
    return mapping_workflow.get_reference_dir_for_name(reference)


### Data input handling independent of wildcards ######################################################################


### Global rule-set stuff #############################################################################################


def get_outputs():
    sample_names = get_sample_names()
    references = get_reference_names()

    outputs = {}
    if "bcftools" in config["variants"]["callers"]:
        step = "filtered" if config["variants__bcftools"]["do_postfilter"] else "all"
        outputs["variants_bcftools"] = expand(
            f"results/variants/{{reference}}/{{sample}}/bcftools_{step}.vcf", sample=sample_names, reference=references
        )
    if "ivar" in config["variants"]["callers"]:
        step = "filtered" if config["variants__ivar"]["do_postfilter"] else "all"
        outputs["variants_ivar"] = expand(
            f"results/variants/{{reference}}/{{sample}}/ivar_{step}.{{ext}}",
            sample=sample_names,
            reference=references,
            ext=["html", "vcf"],
        )
    if "freebayes" in config["variants"]["callers"]:
        step = "filtered" if config["variants__freebayes"]["do_postfilter"] else "all"
        outputs["variants_freebayes"] = expand(
            f"results/variants/{{reference}}/{{sample}}/freebayes_{step}.vcf", sample=sample_names, reference=references
        )
    if "mutect2" in config["variants"]["callers"]:
        step = "filtered" if config["variants__mutect2"]["do_postfilter"] else "all"
        outputs["variants_mutect2"] = expand(
            f"results/variants/{{reference}}/{{sample}}/mutect2_{step}.vcf", sample=sample_names, reference=references
        )

    if "ivar" in config["consensus"]["callers"]:
        outputs["consensus_ivar"] = expand(
            "results/consensus/ivar/{reference}/{sample}.fa", sample=sample_names, reference=references
        )

    if not outputs:
        raise ValueError("No outputs defined for variant calling or consensus calling.")
    return outputs


def infer_reference_dict(wildcards):
    return f"{get_reference_dir_for_name(wildcards.reference)}/{wildcards.reference}.dict"


def infer_reference_fasta(wildcards):
    return f"{get_reference_dir_for_name(wildcards.reference)}/{wildcards.reference}.fa"


def infer_reference_faidx(wildcards):
    return f"{get_reference_dir_for_name(wildcards.reference)}/{wildcards.reference}.fa.fai"


### Parameter parsing from config #####################################################################################


def get_bcftools_mpileup_params():
    extra = [
        "--min-MQ {val}".format(val=config["variants__bcftools"]["min_mapping_quality"]),
        "--min-BQ {val}".format(val=config["variants__bcftools"]["min_base_quality"]),
        "--adjust-MQ {val}".format(val=config["variants__bcftools"]["adjust_mapping_quality"]),
        "--gap-frac {val}".format(val=config["variants__bcftools"]["min_fraction_of_gaps"]),
        "--max-depth {val}".format(val=config["variants__bcftools"]["max_read_depth"]),
        "--max-idepth {val}".format(val=config["variants__bcftools"]["max_depth_indels"]),
        "--max-BQ {val}".format(val=config["variants__bcftools"]["max_base_quality_cap"]),
        "--ext-prob {val}".format(val=config["variants__bcftools"]["gap_ext_prob"]),
        "--open-prob {val}".format(val=config["variants__bcftools"]["open_prob"]),
    ]
    if config["variants__bcftools"]["count_orphans"]:
        extra.append("--count-orphans")
    match config["variants__bcftools"]["compute_base_alignment_quality"]:
        case "no":
            extra.append("--no-BAQ")
        case "redo":
            extra.append("--redo-BAQ")
        case "full":
            extra.append("--full-BAQ")
        case _:
            raise ValueError("Unexpected value for variants__bcftools->compute_base_alignment_quality.")

    if fields := config["variants__bcftools"]["mpileup_annotate"]:
        extra.append(f"--annotate {','.join(fields)}")
    return " ".join(extra)


def get_bcftools_calling_params():
    extra = []
    if config["variants__bcftools"]["keep_alts"]:
        extra.append("--keep-alts")
    if config["variants__bcftools"]["variants_only"]:
        extra.append("--variants-only")
    extra.append(f'--ploidy {config["variants__bcftools"]["ploidy"]}')
    extra.append(f'--prior {config["variants__bcftools"]["prior"]}')
    if fields := config["variants__bcftools"]["caller_annotate"]:
        extra.append(f"--annotate {','.join(fields)}")
    return " ".join(extra)


def get_bcftools_norm_params(tool: str):
    if not config[f"variants__{tool}"]["do_postfilter"]:
        return ""

    cmd = []
    if form := config[f"variants__{tool}"]["postfilter"]["multiallelic"]:
        cmd.append(f" --multiallelics {form}")

    if config[f"variants__{tool}"]["postfilter"]["atomize"]:
        cmd.append(" --atomize")

        if kind := config[f"variants__{tool}"]["postfilter"]["atom_overlaps"]:
            x = "." if kind == "missing" else "\*"
            cmd.append(f" --atom-overlaps {x}")
    if not cmd:
        raise ValueError("No normalization options selected for {tool}.")
    cmd.append("--check-ref w")
    return " ".join(cmd)


def parse_tags_for_bcftools_fill_tags(tool: str):
    if not config[f"variants__{tool}"]["do_postfilter"]:
        return ""
    tags = config[f"variants__{tool}"]["postfilter"]["additional_tags"]
    if tags:
        return f"--tags {tags}"
    else:
        logger.warning(
            "WARNING: No additional tags were specified for {tool}, but this is illogicaly translated by bcftools to all possible tags."
        )
        return ""


def get_bcftools_filter_params(tool: str):
    # common function to parse tool filtering params
    if not config[f"variants__{tool}"]["do_postfilter"]:
        return ""

    extra = []
    if (value := config[f"variants__{tool}"]["postfilter"]["snp_gap"]) is not None:
        extra.append(f"--SnpGap {value}")
    if (value := config[f"variants__{tool}"]["postfilter"]["indel_gap"]) is not None:
        extra.append(f"--IndelGap {value}")
    if (value := config[f"variants__{tool}"]["postfilter"]["set_gts_for_failed"]) is not None:
        x = "'.'" if value == "missing" else "0"
        extra.append(f"--set-GTs {x}")
    if value := config[f"variants__{tool}"]["postfilter"]["include"]:
        extra.append(f"--include '{value}'")
    if value := config[f"variants__{tool}"]["postfilter"]["exclude"]:
        extra.append(f"--exclude '{value}'")
    return " ".join(extra)


def get_bcftools_view_filter_params(tool: str):
    if not config[f"variants__{tool}"]["do_postfilter"]:
        return ""
    return "--trim-alt-alleles" if config[f"variants__{tool}"]["postfilter"]["trim_alt_alleles"] else ""


def parse_samtools_mpileup_for_ivar(kind: str):
    samtools_params = [
        "--max-depth {value}".format(value=config[f"{kind}__ivar"]["max_read_depth"]),
        "--min-MQ {value}".format(value=config[f"{kind}__ivar"]["min_mapping_quality"]),
        "--min-BQ {value}".format(value=config[f"{kind}__ivar"]["min_base_quality"]),
    ]

    if config[f"{kind}__ivar"]["count_orphans"]:
        samtools_params.append("--count-orphans")
    if config[f"{kind}__ivar"]["absolutely_all_positions"]:
        samtools_params.append("-aa")
    if config[f"{kind}__ivar"]["redo_base_alignment_quality"]:
        samtools_params.append("--redo-BAQ")
    else:
        samtools_params.append("--no-BAQ")

    return " ".join(samtools_params)


def parse_ivar_params_for_variants():
    return " ".join(
        [
            "-q {value}".format(value=config["variants__ivar"]["min_base_quality_threshold"]),
            "-t {value}".format(value=config["variants__ivar"]["min_frequency_threshold"]),
            "-m {value}".format(value=config["variants__ivar"]["min_read_depth"]),
        ]
    )


def parse_ivar_params_for_consensus():
    return " ".join(
        [
            "-q {value}".format(value=config["consensus__ivar"]["consensus_base_quality_threshold"]),
            "-t {value}".format(value=config["consensus__ivar"]["consensus_frequency_threshold"]),
            "-m {value}".format(value=config["consensus__ivar"]["min_consensus_depth"]),
            "-c {value}".format(value=config["consensus__ivar"]["insertion_frequency_threshold"]),
        ]
    )


def parse_freebayes_params():
    extra = []

    if val := config["variants__freebayes"]["min_base_quality"]:
        extra.append(f"--min-base-quality {val}")

    if val := config["variants__freebayes"]["min_mapping_quality"]:
        extra.append(f"--min-mapping-quality {val}")

    if val := config["variants__freebayes"]["mismatch_base_quality_threshold"]:
        extra.append(f"--mismatch-base-quality-threshold {val}")

    if val := config["variants__freebayes"]["read_max_mismatch_fraction"]:
        extra.append(f"--read-max-mismatch-fraction {val}")

    if (val := config["variants__freebayes"]["read_mismatch_limit"]) is not None:
        extra.append(f"--read-mismatch-limit {val}")

    if (val := config["variants__freebayes"]["read_indel_limit"]) is not None:
        extra.append(f"--read-indel-limit {val}")

    if (val := config["variants__freebayes"]["read_snp_limit"]) is not None:
        extra.append(f"--read-snp-limit {val}")

    if val := config["variants__freebayes"]["indel_exclusion_window"]:
        extra.append(f"--indel-exclusion-window {val}")

    if val := config["variants__freebayes"]["min_alternate_fraction"]:
        extra.append(f"--min-alternate-fraction {val}")

    if val := config["variants__freebayes"]["min_alternate_count"]:
        extra.append(f"--min-alternate-count {val}")

    if val := config["variants__freebayes"]["min_coverage"]:
        extra.append(f"--min-coverage {val}")

    if (val := config["variants__freebayes"]["limit_coverage"]) is not None:
        extra.append(f"--limit-coverage {val}")

    if (val := config["variants__freebayes"]["skip_coverage"]) is not None:
        extra.append(f"--skip-coverage {val}")

    if val := config["variants__freebayes"]["ploidy"]:
        extra.append(f"--ploidy {val}")

    if config["variants__freebayes"]["gvcf"]:
        extra.append("--gvcf")

    if config["variants__freebayes"]["report_all_haplotype_alleles"]:
        extra.append("--report-all-haplotype-alleles")

    if config["variants__freebayes"]["report_monorphic"]:
        extra.append("--report-monomorphic")

    return " ".join(extra)


def parse_mutect2_call_params():
    extra = [
        "--create-output-variant-index",
        f'--base-quality-score-threshold {config["variants__mutect2"]["base_quality_score_threshold"]}',
        f'--callable-depth {config["variants__mutect2"]["callable_depth"]}',
        f'--f1r2-max-depth {config["variants__mutect2"]["f1r2_max_depth"]}',
        f'--f1r2-median-mq {config["variants__mutect2"]["f1r2_median_mq"]}',
        f'--f1r2-min-bq {config["variants__mutect2"]["f1r2_min_bq"]}',
        f'--min-base-quality-score {config["variants__mutect2"]["min_base_quality_score"]}',
        f'--max-reads-per-alignment-start {config["variants__mutect2"]["max_reads_per_alignment_start"]}',
    ]
    return " ".join(extra)


def parse_mutect2_filter_params():
    if not config["variants__mutect2"]["do_postfilter"]:
        return ""

    cfg = config["variants__mutect2"]["postfilter"]
    extra = [
        "--create-output-variant-index",
        f'--min-median-base-quality {cfg["min_median_base_quality"]}',
        f'--min-median-mapping-quality {cfg["min_median_mapping_quality"]}',
        f'--min-median-read-position {cfg["min_median_read_position"]}',
        f'--min-reads-per-strand {cfg["min_reads_per_strand"]}',
        f'--min-slippage-length {cfg["min_slippage_length"]}',
        f'--pcr-slippage-rate {cfg["pcr_slippage_rate"]}',
    ]
    if cfg["microbial_mode"]:
        extra.append("--microbial-mode")
    return " ".join(extra)


def get_all_relevant_extra_params():
    extra = ""

    if "freebayes" in config["variants"]["callers"]:
        extra += f"Freebayes: {parse_freebayes_params()}\n"
    if "mutect2" in config["variants"]["callers"]:
        extra += f"Mutect2 call: {parse_mutect2_call_params()}\n"
        extra += f"Mutect2 filter: {parse_mutect2_filter_params()}\n"
    if "ivar" in config["variants"]["callers"]:
        extra += f"IVAR variants:\n"
        extra += f"\tsamtools mpileup: {parse_samtools_mpileup_for_ivar('variants')}\n"
        extra += f"\tIVAR variants: {parse_ivar_params_for_variants()}\n"
    if "bcftools" in config["variants"]["callers"]:
        extra += f"BCFtools variants:\n"
        extra += f"\tmpileup: {get_bcftools_mpileup_params()}\n"
        extra += f"\tcall: {get_bcftools_calling_params()}\n"
        extra += f"\tnorm: {get_bcftools_norm_params('bcftools')}\n"
        extra += f"\tfilter: {get_bcftools_filter_params('bcftools')}\n"
        extra += f"\tview filter: {get_bcftools_view_filter_params('bcftools')}\n"
    if "ivar" in config["consensus"]["callers"]:
        extra += "IVAR consensus:\n"
        extra += f"\tsamtools mpileup: {parse_samtools_mpileup_for_ivar('consensus')}\n"
        extra += f"\tIVAR consensus: {parse_ivar_params_for_consensus()}\n"
    return extra


### Resource handling #################################################################################################


def get_threads_for_bcftools():
    return min(config["threads"]["variants__bcftools"], config["max_threads"])


def get_threads_for_freebayes():
    return min(config["threads"]["variants__freebayes"], config["max_threads"])


def get_threads_for_mutect2():
    return min(config["threads"]["variants__mutect2"], config["max_threads"])


def get_mem_mb_for_freebayes(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["variants__freebayes_mem_mb"] * attempt)


def get_mem_mb_for_mutect2(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["variants__mutect2_mem_mb"] * attempt)
