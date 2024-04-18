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
        if config["variants__bcftools"]["do_postfilter"]:
            outputs["variants_bcftools"] = expand(
                "results/variants_bcftools/{reference}/filtered/{sample}.vcf", sample=sample_names, reference=references
            )
        else:
            outputs["variants_bcftools"] = expand(
                "results/variants_bcftools/{reference}/normalized/{sample}.vcf",
                sample=sample_names,
                reference=references,
            )
    if "ivar" in config["variants"]["callers"]:
        if config["variants__ivar"]["do_postfilter"]:
            outputs["variants_ivar"] = expand(
                "results/variants_ivar/{reference}/{sample}/filtered.html", sample=sample_names, reference=references
            )
        else:
            outputs["variants_ivar"] = expand(
                "results/variants_ivar/{reference}/{sample}/all.vcf", sample=sample_names, reference=references
            )
    if "freebayes" in config["variants"]["callers"]:
        if config["variants__freebayes"]["do_postfilter"]:
            outputs["variants_freebayes"] = expand(
                "results/variants_freebayes/{reference}/{sample}/filtered.vcf",
                sample=sample_names,
                reference=references,
            )
        else:
            outputs["variants_freebayes"] = expand(
                "results/variants_freebayes/{reference}/{sample}/original.vcf",
                sample=sample_names,
                reference=references,
            )

    if "mutect2" in config["variants"]["callers"]:
        outputs["variants_mutect2"] = expand(
            "results/variants_mutect2/{reference}/{sample}/filtered.vcf",
            sample=sample_names,
            reference=references,
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
        "--adjust-MQ {val}".format(val=config["variants__bcftools"]["adjust_quality_coef"]),
        "--gap-frac {val}".format(val=config["variants__bcftools"]["min_fraction_of_gaps"]),
    ]
    if config["variants__bcftools"]["redo_base_alignment_qual"]:
        extra.append("--redo-BAQ")
    if form := config["variants__bcftools"]["annotation"]:
        extra.append(f"--annotate {form}")
    return " ".join(extra)


def get_bcftools_calling_params():
    extra = []
    if config["variants__bcftools"]["keep_alts"]:
        extra.append("--keep-alts")
    if config["variants__bcftools"]["variants_only"]:
        extra.append("--variants-only")
    if form := config["variants__bcftools"]["ploidy"]:
        extra.append(f"--ploidy {form}")
    if fields := config["variants__bcftools"]["format_fields"]:
        extra.append(f"--format-fields {fields}")
    return " ".join(extra)


def get_bcftools_norm_params(tool: str):
    # common function to parse tool normalization params
    if not config[f"variants__{tool}"]["do_postfilter"]:
        return ""

    cmd = "--check-ref w"
    if form := config[f"variants__{tool}"]["postfilter"]["multiallelic"]:
        cmd += f" --multiallelics {form}"

        if form := config[f"variants__{tool}"]["postfilter"]["multi_overlaps"]:
            x = "'.'" if form == "missing" else "0"
            cmd += f" --multi-overlaps {x}"

    if config[f"variants__{tool}"]["postfilter"]["atomize"]:
        x = "'.'" if form == "missing" else "0"
        cmd += f" --multi-overlaps {x}"

        if kind := config[f"variants__{tool}"]["postfilter"]["atom_overlaps"]:
            x = "." if kind == "missing" else "\*"
            cmd += f" --atom-overlaps {x}"
    return cmd


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


def parse_samtools_params_for_variants():
    samtools_params = ["-aa --no-BAQ"]
    if config["variants__ivar"]["count_orphans"]:
        samtools_params.append("--count-orphans")

    samtools_params.append("--max-depth {value}".format(value=config["variants__ivar"]["max_read_depth"]))
    samtools_params.append("--min-MQ {value}".format(value=config["variants__ivar"]["min_mapping_quality"]))
    samtools_params.append("--min-BQ {value}".format(value=config["variants__ivar"]["min_base_quality"]))
    return " ".join(samtools_params)


def parse_ivar_params_for_variants():
    ivar_params = []

    ivar_params.append("-q {value}".format(value=config["variants__ivar"]["min_base_quality_threshold"]))
    ivar_params.append("-t {value}".format(value=config["variants__ivar"]["min_frequency_threshold"]))
    ivar_params.append("-m {value}".format(value=config["variants__ivar"]["min_read_depth"]))
    return " ".join(ivar_params)


def parse_samtools_params_for_consensus():
    samtools_params = ["--no-BAQ"]
    if config["consensus__ivar"]["count_orphans"]:
        samtools_params.append("--count-orphans")

    samtools_params.append("--max-depth {value}".format(value=config["consensus__ivar"]["max_read_depth"]))
    samtools_params.append("--min-MQ {value}".format(value=config["consensus__ivar"]["min_mapping_quality"]))
    samtools_params.append("--min-BQ {value}".format(value=config["consensus__ivar"]["min_base_quality"]))
    return " ".join(samtools_params)


def parse_ivar_params_for_consensus():
    ivar_params = []

    ivar_params.append("-q {value}".format(value=config["consensus__ivar"]["consensus_base_quality_threshold"]))
    ivar_params.append("-t {value}".format(value=config["consensus__ivar"]["consensus_frequency_threshold"]))
    ivar_params.append("-m {value}".format(value=config["consensus__ivar"]["min_consensus_depth"]))

    return " ".join(ivar_params)


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


### Resource handling #################################################################################################


def get_threads_for_bcftools():
    return min(config["threads"]["variants__bcftools"], config["max_threads"])


def get_threads_for_freebayes():
    return min(config["threads"]["variants__freebayes"], config["max_threads"])


def get_mem_mb_for_freebayes(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["variants__freebayes_mem_mb"] * attempt)
