import os


include: "constants.smk"


rule create_bedfile_asm:
    input:
        fa_asm=os.path.join(INPUT_DIR_ASM, "{sm}.fa.gz"),
    output:
        fai=os.path.join(INPUT_DIR_ASM, "{sm}.fa.gz.fai"),
        bedfile=os.path.join(OUTPUT_DIR, "calculate_cov", "asm_{sm}_wg.bed"),
        json=os.path.join(OUTPUT_DIR, "calculate_cov", "asm_{sm}_wg.json"),
    conda:
        "../envs/tools.yaml"
    resources:
        mem=4,
    log:
        os.path.join(LOGS_DIR, "calculate_cov", "create_bedfile_{sm}.log"),
    benchmark:
        os.path.join(BENCHMARKS_DIR, "calculate_cov", "create_bedfile_{sm}.tsv")
    shell:
        """
        # create fasta index for the input diploid assembly
        samtools faidx {input.fa_asm} 2> {log}

        # create a bed file for whole genome
        {{ cat {input.fa_asm}.fai | awk '{{ print $1"\\t0\\t"$2 }}' | bedtools sort -i - ;}} > {output.bedfile} 2>> {log}

        # make a json file pointing to the asm_wg.bed
        echo "{{" > {output.json}
        echo \"asm_wg\" : \"asm_wg.bed\" >> {output.json}
        echo "}}" >> {output.json}
        """


rule bam2cov:
    input:
        bamfile=os.path.join(OUTPUT_DIR, "secphase", "{sm}", "{sm}_corrected.bam"),
        bedfile_json=rules.create_bedfile_asm.output.json,
    output:
        os.path.join(OUTPUT_DIR, "calculate_cov", "read_alignment_{sm}.cov"),
    params:
        output_type="c",
        output_format="only_total",
    singularity:
        "docker://mobinasri/flagger:v0.4.0"
    resources:
        mem=config["calculate_cov"]["mem"],
    threads: config["calculate_cov"]["threads"]
    log:
        os.path.join(LOGS_DIR, "calculate_cov", "bam2cov_{sm}.log"),
    benchmark:
        os.path.join(BENCHMARKS_DIR, "calculate_cov", "bam2cov_{sm}.tsv")
    shell:
        """
        bam2cov \
        -i {input.bamfile} \
        -O "{params.output_type}" \
        -o {output} \
        -j {input.bedfile_json} \
        -t 8 \
        -f {params.output_format} 2> {log}
        """


rule calculate_cov_all:
    input:
        expand(rules.create_bedfile_asm.output, sm=SAMPLE_NAMES),
        expand(rules.bam2cov.output, sm=SAMPLE_NAMES),
