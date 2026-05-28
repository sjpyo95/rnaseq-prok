rule fastqc_raw:
    input:
        fq1 = get_fq1,
        fq2 = get_fq2,
    output:
        html_r1 = "results/qc/raw/{sample}_R1_fastqc.html",
        html_r2 = "results/qc/raw/{sample}_R2_fastqc.html",
        zip_r1  = "results/qc/raw/{sample}_R1_fastqc.zip",
        zip_r2  = "results/qc/raw/{sample}_R2_fastqc.zip",
    params:
        outdir = "results/qc/raw",
    threads: config["threads"]["fastqc"]
    conda: "../../envs/qc.yaml"
    log: "logs/fastqc_raw/{sample}.log"
    shell:
        "fastqc -t {threads} {input.fq1} {input.fq2} -o {params.outdir} > {log} 2>&1"


rule fastqc_trimmed:
    input:
        fq1 = "results/trim/{sample}_R1.fastq.gz",
        fq2 = "results/trim/{sample}_R2.fastq.gz",
    output:
        html_r1 = "results/qc/trimmed/{sample}_R1_fastqc.html",
        html_r2 = "results/qc/trimmed/{sample}_R2_fastqc.html",
        zip_r1  = "results/qc/trimmed/{sample}_R1_fastqc.zip",
        zip_r2  = "results/qc/trimmed/{sample}_R2_fastqc.zip",
    params:
        outdir = "results/qc/trimmed",
    threads: config["threads"]["fastqc"]
    conda: "../../envs/qc.yaml"
    log: "logs/fastqc_trimmed/{sample}.log"
    shell:
        "fastqc -t {threads} {input.fq1} {input.fq2} -o {params.outdir} > {log} 2>&1"


rule multiqc:
    input:
        expand("results/qc/raw/{sample}_R1_fastqc.zip",     sample=SAMPLES),
        expand("results/qc/raw/{sample}_R2_fastqc.zip",     sample=SAMPLES),
        expand("results/qc/trimmed/{sample}_R1_fastqc.zip", sample=SAMPLES),
        expand("results/qc/trimmed/{sample}_R2_fastqc.zip", sample=SAMPLES),
        expand("results/trim/{sample}_fastp.json",          sample=SAMPLES),
        expand("logs/align/{sample}.log",                   sample=SAMPLES),
    output:
        "results/qc/multiqc_report.html",
    params:
        outdir = "results/qc",
        dirs   = "results/qc results/trim logs/align",
    conda: "../../envs/qc.yaml"
    log: "logs/multiqc.log"
    shell:
        "multiqc {params.dirs} -o {params.outdir} --force > {log} 2>&1"
