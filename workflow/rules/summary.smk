rule concat__consensuses_for_references:
    input:
        consensuses=infer_consensuses_for_reference_tool,
    output:
        report(
            "results/_aggregation/consensus/{reference}_{tool}.fa",
            category="_Summary",
            labels={
                "Reference": "{reference}",
                "Type": "Consensus - {tool}",
            },
        ),
    log:
        "logs/concat/consensuses_for_references/{reference}_{tool}.log",
    conda:
        "../envs/awk_sed.yaml"
    localrule: True
    shell:
        "cat {input.consensuses} 1> {output} 2> {log}"
