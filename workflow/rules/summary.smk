rule concat__consensuses:
    input:
        expand("results/consensus/{{reference}}/{sample}/ivar.fa", sample=get_sample_names()),
    output:
        report(
            "results/_aggregation/consensus/{reference}.fa",
            category="_aggregation",
            labels={
                "Type": "Consensuses for {reference}",
            },
        ),
    log:
        "logs/concat/consensuses_for_references/{reference}.log",
    conda:
        "../envs/awk_sed.yaml"
    localrule: True
    shell:
        "cat {input} 1> {output} 2> {log}"
