rule fastqc_raw:
    input:
        fq1 = get_fq1,
        fq2 = get_fq2,
    output:
        html_r1 = OUTDIR + "/qc/raw/{sample}_R1_fastqc.html",
        html_r2 = OUTDIR + "/qc/raw/{sample}_R2_fastqc.html",
        zip_r1  = OUTDIR + "/qc/raw/{sample}_R1_fastqc.zip",
        zip_r2  = OUTDIR + "/qc/raw/{sample}_R2_fastqc.zip",
    params:
        outdir = OUTDIR + "/qc/raw",
    threads: config["threads"]["fastqc"]
    conda: "../../envs/qc.yaml"
    log: "logs/fastqc_raw/{sample}.log"
    shell:
        "fastqc -t {threads} {input.fq1} {input.fq2} -o {params.outdir} > {log} 2>&1"


rule fastqc_trimmed:
    input:
        fq1 = OUTDIR + "/trim/{sample}_R1.fastq.gz",
        fq2 = OUTDIR + "/trim/{sample}_R2.fastq.gz",
    output:
        html_r1 = OUTDIR + "/qc/trimmed/{sample}_R1_fastqc.html",
        html_r2 = OUTDIR + "/qc/trimmed/{sample}_R2_fastqc.html",
        zip_r1  = OUTDIR + "/qc/trimmed/{sample}_R1_fastqc.zip",
        zip_r2  = OUTDIR + "/qc/trimmed/{sample}_R2_fastqc.zip",
    params:
        outdir = OUTDIR + "/qc/trimmed",
    threads: config["threads"]["fastqc"]
    conda: "../../envs/qc.yaml"
    log: "logs/fastqc_trimmed/{sample}.log"
    shell:
        "fastqc -t {threads} {input.fq1} {input.fq2} -o {params.outdir} > {log} 2>&1"


rule multiqc:
    input:
        expand(OUTDIR + "/qc/raw/{sample}_R1_fastqc.zip",     sample=SAMPLES),
        expand(OUTDIR + "/qc/raw/{sample}_R2_fastqc.zip",     sample=SAMPLES),
        expand(OUTDIR + "/qc/trimmed/{sample}_R1_fastqc.zip", sample=SAMPLES),
        expand(OUTDIR + "/qc/trimmed/{sample}_R2_fastqc.zip", sample=SAMPLES),
        expand(OUTDIR + "/trim/{sample}_fastp.json",           sample=SAMPLES),
        expand("logs/align/{sample}.log",                      sample=SAMPLES),
    output:
        OUTDIR + "/qc/multiqc_report.html",
    params:
        outdir = OUTDIR + "/qc",
        dirs   = OUTDIR + "/qc " + OUTDIR + "/trim logs/align",
    conda: "../../envs/qc.yaml"
    log: "logs/multiqc.log"
    shell:
        "multiqc {params.dirs} -o {params.outdir} --force > {log} 2>&1"
