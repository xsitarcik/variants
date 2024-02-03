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


def get_bcftools_norm_params():
    if not config["variants__bcftools"]["do_postfilter"]:
        return ""

    base = "--check-ref w"
    if form := config["variants__bcftools"]["postfilter"]["multiallelic"]:
        return base + f" --multiallelics {form}"
    return base


def get_bcftools_filter_params():
    if not config["variants__bcftools"]["do_postfilter"]:
        return ""

    extra = []
    if value := config["variants__bcftools"]["postfilter"]["snp_gap"]:
        extra.append(f"--SnpGap {value}")
    if value := config["variants__bcftools"]["postfilter"]["indel_gap"]:
        extra.append(f"--IndelGap {value}")
    if value := config["variants__bcftools"]["postfilter"].get("set_gts_for_failed", False):
        extra.append(f"--set-GTs {value}")
    if value := config["variants__bcftools"]["postfilter"].get("include", False):
        extra.append(f"--include '{value}'")
    if value := config["variants__bcftools"]["postfilter"].get("exclude", False):
        extra.append(f"--exclude '{value}'")
    return " ".join(extra)


def get_bcftools_view_filter_params():
    if not config["variants__bcftools"]["do_postfilter"]:
        return ""
    return "--trim-alt-alleles" if config["variants__bcftools"]["postfilter"]["trim_alt_alleles"] else ""


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


### Resource handling #################################################################################################


def get_threads_for_bcftools():
    return min(config["threads"]["variants__bcftools"], config["max_threads"])
