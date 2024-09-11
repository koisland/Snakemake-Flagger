import os


include: "constants.smk"


rule flagger:
    input:
        cov=os.path.join(OUTPUT_DIR, "calculate_cov", "read_alignment_{sm}.cov"),
        alpha=config["flagger"]["alpha"],
    output:
        directory(os.path.join(OUTPUT_DIR, "flagger", "{sm}")),
    singularity:
        "docker://mobinasri/flagger:v1.0.0"
    resources:
        mem=config["flagger"]["mem"],
    threads: 16
    params:
        label_names="Err,Dup,Hap,Col",
        track_name="hmm_flagger_v1.0",
        num_iter=100,
        chunk_len=20000000,
        window_len=4000,
        convergence_tolerance=0.001,
        high_mapq_ratio=0.25,
        model_type="gaussian",
    log:
        os.path.join(LOGS_DIR, "flagger", "flagger_{sm}.log"),
    benchmark:
        os.path.join(BENCHMARKS_DIR, "flagger", "flagger_{sm}.tsv")
    shell:
        """
        mkdir -p {output}
        hmm_flagger \
        --input {input.cov} \
        --alpha {input.alpha} \
        --outputDir {output}  \
        --labelNames {params.label_names} \
        --chunkLen {params.chunk_len} \
        --windowLen {params.window_len} \
        --maxHighMapqRatio {params.high_mapq_ratio} \
        --modelType {params.model_type} \
        --iterations {params.num_iter} \
        --trackName {params.track_name} \
        --convergenceTol {params.convergence_tolerance} \
        --threads {threads} 2> {log}
        """


rule flagger_all:
    input:
        expand(rules.flagger.output, sm=SAMPLE_NAMES),
