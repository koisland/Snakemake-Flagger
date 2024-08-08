import os


include: "constants.smk"


checkpoint cov2counts:
    input:
        os.path.join(OUTPUT_DIR, "calculate_cov", "read_alignment_{sm}.cov"),
    output:
        os.path.join(OUTPUT_DIR, "cov_dist_fit", "read_alignment_{sm}.counts"),
    singularity:
        "docker://mobinasri/flagger:v0.4.0"
    resources:
        mem=config["cov_dist_fit"]["mem"],
    log:
        os.path.join(LOGS_DIR, "cov_dist_fit", "cov2counts_{sm}.log"),
    benchmark:
        os.path.join(BENCHMARKS_DIR, "cov_dist_fit", "cov2counts_{sm}.tsv")
    shell:
        """
        cov2counts -i {input} -o {output} 2> {log}
        """


# TODO: Need to get {EXPECTED_COVERAGE}.
def get_expected_coverage(wc):
    raise NotImplementedError()

    read_counts = checkpoints.cov2counts.get(**wc).output
    with open(read_counts, "rt") as fh:
        pass


rule fit_gmm:
    input:
        rules.cov2counts.output,
    output:
        os.path.join(OUTPUT_DIR, "cov_dist_fit", "read_alignment_{sm}.table"),
    singularity:
        "docker://mobinasri/flagger:v0.4.0"
    resources:
        mem=config["cov_dist_fit"]["mem"],
    log:
        os.path.join(LOGS_DIR, "cov_dist_fit", "fit_gmm_{sm}.log"),
    benchmark:
        os.path.join(BENCHMARKS_DIR, "cov_dist_fit", "fit_gmm_{sm}.tsv")
    params:
        expected_coverage=0.2,
    shell:
        """
        python3 /home/programs/src/fit_gmm.py \
        --counts {input} \
        --cov {params.expected_coverage} \
        --output {output} 2> {log}
        """


rule cov_dist_fit_all:
    input:
        expand(rules.cov2counts.output, sm=SAMPLE_NAMES),
        expand(rules.fit_gmm.output, sm=SAMPLE_NAMES),
