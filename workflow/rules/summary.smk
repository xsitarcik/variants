rule concat__consensuses_for_references:
    input:
        consensuses=infer_consensuses_for_reference,
    output:
        report(
            "results/_aggregation/consensus/{reference}.fa",
            category="Summary",
            labels={
                "Type": "Consensus for {reference}",
            },
        ),
    log:
        "logs/concat/consensuses_for_references/{reference}.log",
    conda:
        "../envs/awk_sed.yaml"
    localrule: True
    shell:
        "cat {input.consensuses} 1> {output} 2> {log}"
