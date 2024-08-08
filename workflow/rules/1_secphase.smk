import os


include: "constants.smk"


rule secphase:
    input:
        fa_asm=os.path.join(INPUT_DIR_ASM, "{sm}.fa.gz"),
        # Secondary reads should be included
        bam_sorted=os.path.join(INPUT_DIR_BAM, "{sm}.bam"),
    output:
        os.path.join(OUTPUT_DIR, "secphase", "{sm}", "{sm}.out.log"),
    singularity:
        "docker://mobinasri/secphase:v0.4.3"
    threads: config["secphase"]["threads"]
    resources:
        mem=config["secphase"]["mem"],
    log:
        os.path.join(LOGS_DIR, "secphase", "secphase_{sm}.log"),
    benchmark:
        os.path.join(BENCHMARKS_DIR, "secphase", "secphase_{sm}.tsv")
    shell:
        """
        secphase --hifi \
        -i {input.bam_sorted} \
        -f {input.fa_asm} \
        --outDir {output} \
        --prefix {wildcards.sm} \
        --threads {threads} 2> {log}
        """


rule correct_bam:
    input:
        bam_sorted=rules.secphase.input.bam_sorted,
        phase_log=rules.secphase.output,
    output:
        corrected_bam=os.path.join(OUTPUT_DIR, "secphase", "{sm}", "{sm}_corrected.bam"),
    params:
        primary_only="--primaryOnly",
    singularity:
        "docker://mobinasri/secphase:v0.4.3"
    resources:
        mem=config["secphase"]["mem"],
    log:
        os.path.join(LOGS_DIR, "secphase", "correct_bam_{sm}.log"),
    benchmark:
        os.path.join(BENCHMARKS_DIR, "secphase", "correct_bam_{sm}.tsv")
    shell:
        """
        correct_bam \
        -i {input.bam_sorted} \
        -P {input.phase_log} \
        -o {output.corrected_bam} \
        {params.primary_only} 2> {log}
        """


rule secphase_all:
    input:
        expand(rules.secphase.output, sm=SAMPLE_NAMES),
        expand(rules.correct_bam.output, sm=SAMPLE_NAMES),
