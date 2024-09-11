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


# Remove N's
rule create_wg_bedfile:
    input:
        fa_asm=os.path.join(INPUT_DIR_ASM, "{sm}.fa.gz"),
        fai=rules.index_asm.output,
    output:
        fa_asm=temp(os.path.join(INPUT_DIR_ASM, "{sm}.fa")),
        wg_bed=os.path.join(OUTPUT_DIR, "calculate_cov", "asm_{sm}_wg.bed"),
    singularity:
        "docker://mobinasri/flagger:v1.0.0"
    resources:
        mem=4,
    log:
        os.path.join(LOGS_DIR, "calculate_cov", "create_wg_bedfile_{sm}.log"),
    shell:
        """
        zcat {input.fa_asm} > {output.fa_asm}
        {{ awk -v OFS="\\t" '{{ print $1, 0, $2}}' {input.fai} | \
            bedtools sort -i - | \
            bedtools subtract -a - -b <(python3 /home/scripts/get_N_coords.py --inputFasta {output.fa_asm});}} > {output.wg_bed} 2> {log}
        """


rule create_wg_json:
    input:
        rules.create_wg_bedfile.output.wg_bed,
    output:
        json=os.path.join(OUTPUT_DIR, "calculate_cov", "asm_{sm}_wg.json"),
    singularity:
        "docker://mobinasri/flagger:v1.0.0"
    resources:
        mem=4,
    log:
        os.path.join(LOGS_DIR, "calculate_cov", "create_bedfile_{sm}.log"),
    shell:
        """
        # make a json file pointing to the asm_wg.bed
        echo "{{" > {output.json}
        echo \\"asm_wg\\": \\"{input}\\" >> {output.json}
        echo "}}" >> {output.json}
        """


rule bam2cov:
    input:
        bamfile=os.path.join(INPUT_DIR_BAM, "{sm}.bam"),
        bedfile_json=rules.create_wg_json.output.json,
        bias_regions_json=lambda wc: (
            config["calculate_cov"]["biased_regions"][str(wc.sm)]
            if config["calculate_cov"].get("biased_regions", {}).get(str(wc.sm))
            else []
        ),
    output:
        os.path.join(OUTPUT_DIR, "calculate_cov", "read_alignment_{sm}.cov"),
    params:
        baseline_annotation="whole_genome",
        mapq_threshold=20,
        clip_ratio_threshold=0.1,
        downsample_rate=1.0,
        bias_regions_json=lambda wc, input: (
            f"--runBiasDetection --restrictBiasAnnotations {input.bias_regions_json}"
            if input.bias_regions_json
            else ""
        ),
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
        --mapqThreshold {params.mapq_threshold} \
        --clipRatioThreshold {params.clip_ratio_threshold} \
        --downsampleRate {params.downsample_rate} \
        --threads 8 \
        {params.bias_regions_json} \
        --baselineAnnotation {params.baseline_annotation} 2> {log}
        """


rule calculate_cov_all:
    input:
        expand(rules.bam2cov.output, sm=SAMPLE_NAMES),
