import os


include: "constants.smk"


rule flagger:
    input:
        os.path.join(OUTPUT_DIR, "calculate_cov", "read_alignment_{sm}.cov"),
    output:
        directory(os.path.join(OUTPUT_DIR, "flagger", "{sm}")),
    singularity:
        "docker://mobinasri/flagger:v1.0.0"
    resources:
        mem=config["flagger"]["mem"],
    threads:
        16
    params:
        label_names="Err,Dup,Hap,Col",
    log:
        os.path.join(LOGS_DIR, "flagger", "flagger_{sm}.log"),
    benchmark:
        os.path.join(BENCHMARKS_DIR, "flagger", "flagger_{sm}.tsv")
    shell:
        """
        mkdir -p {output}
        hmm_flagger \
        --input {input} \
        --outputDir {output}  \
        --labelNames {params.label_names} \
        --threads {threads} 2> {log}
        """


rule flagger_all:
    input:
        expand(rules.flagger.output, sm=SAMPLE_NAMES),
