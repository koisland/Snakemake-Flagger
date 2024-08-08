import os


include: "constants.smk"


rule find_blocks_from_table:
    input:
        cov=os.path.join(OUTPUT_DIR, "calculate_cov", "read_alignment_{sm}.cov"),
        gmm_fit_table=os.path.join(
            OUTPUT_DIR, "cov_dist_fit", "read_alignment_{sm}.table"
        ),
    output:
        bed_err=os.path.join(OUTPUT_DIR, "extract_blocks", "{sm}_error.bed"),
        bed_dupe=os.path.join(OUTPUT_DIR, "extract_blocks", "{sm}_duplicated.bed"),
        bed_hap=os.path.join(OUTPUT_DIR, "extract_blocks", "{sm}_haploid.bed"),
        bed_collapse=os.path.join(OUTPUT_DIR, "extract_blocks", "{sm}_collapsed.bed"),
    resources:
        mem=config["extract_blocks"]["mem"],
    singularity:
        "docker://mobinasri/flagger:v0.4.0"
    params:
        output_prefix=lambda wc, output: os.path.join(
            os.path.dirname(output.bed_err), f"{wc.sm}"
        ),
    log:
        os.path.join(LOGS_DIR, "extract_blocks", "find_blocks_from_table_{sm}.log"),
    benchmark:
        os.path.join(
            BENCHMARKS_DIR, "extract_blocks", "find_blocks_from_table_{sm}.tsv"
        )
    shell:
        """
        find_blocks_from_table \
        -c {input.cov} \
        -t {input.gmm_fit_table} \
        -p {params.output_prefix} 2> {log}
        """


rule extract_blocks_all:
    input:
        expand(rules.find_blocks_from_table.output, sm=SAMPLE_NAMES),
