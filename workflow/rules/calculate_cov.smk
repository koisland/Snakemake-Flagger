import os


include: "constants.smk"


rule index_asm:
    input:
        fa_asm=os.path.join(INPUT_DIR_ASM, "{sm}.fa.gz"),
    output:
        fai=os.path.join(INPUT_DIR_ASM, "{sm}.fa.gz.fai"),
    singularity:
        "docker://mobinasri/flagger:v1.0.0"
    resources:
        mem=4,
    log:
        os.path.join(LOGS_DIR, "calculate_cov", "index_asm_{sm}.log"),
    shell:
        """
        # create fasta index for the input diploid assembly
        samtools faidx {input.fa_asm} 2> {log}
        """


rule create_bedfile_asm:
    input:
        rules.index_asm.output
    output:
        bedfile=os.path.join(OUTPUT_DIR, "calculate_cov", "asm_{sm}_wg.bed"),
        json=os.path.join(OUTPUT_DIR, "calculate_cov", "asm_{sm}_wg.json"),
    singularity:
        "docker://mobinasri/flagger:v1.0.0"
    resources:
        mem=4,
    log:
        os.path.join(LOGS_DIR, "calculate_cov", "create_bedfile_{sm}.log"),
    shell:
        """
        # create a bed file for whole genome
        {{ cat {input} | awk '{{ print $1"\\t0\\t"$2 }}' | bedtools sort -i - ;}} > {output.bedfile} 2>> {log}

        # make a json file pointing to the asm_wg.bed
        echo "{{" > {output.json}
        echo \\"asm_wg\\": \\"{output.bedfile}\\" >> {output.json}
        echo "}}" >> {output.json}
        """


rule bam2cov:
    input:
        bamfile=os.path.join(INPUT_DIR_BAM, "{sm}.bam"),
        bedfile_json=rules.create_bedfile_asm.output.json,
        bias_regions_json=lambda wc: (
            config["calculate_cov"]["biased_regions"][str(wc.sm)]
            if config["calculate_cov"].get("biased_regions", {}).get(str(wc.sm))
            else []
        )
    output:
        os.path.join(OUTPUT_DIR, "calculate_cov", "read_alignment_{sm}.cov"),
    params:
        baseline_annotation="whole_genome",
        bias_regions_json=lambda wc, input: (
            f"--runBiasDetection {input.bias_regions_json}"
            if input.bias_regions_json else ""
        )
    singularity:
        "docker://mobinasri/flagger:v1.0.0"
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
        --bam {input.bamfile} \
        --output {output} \
        --annotationJson {input.bedfile_json} \
        --threads 8 \
        {params.bias_regions_json} \
        --baselineAnnotation {params.baseline_annotation} 2> {log}
        """


rule calculate_cov_all:
    input:
        expand(rules.create_bedfile_asm.output, sm=SAMPLE_NAMES),
        expand(rules.bam2cov.output, sm=SAMPLE_NAMES),
