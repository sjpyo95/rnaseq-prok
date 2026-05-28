_STRAND_TO_HISAT2 = {"FR": "FR", "RF": "RF", "unstranded": ""}


def get_hisat2_strand_flag(wildcards):
    strand = open(f"results/strand/{wildcards.sample}.strand").read().strip()
    mapped = _STRAND_TO_HISAT2[strand]
    return f"--rna-strandness {mapped}" if mapped else ""


rule hisat2_build:
    """Build HISAT2 index from reference genome."""
    input:
        config["genome"],
    output:
        multiext(config["hisat2_index"],
                 ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",
                 ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"),
    params:
        prefix = config["hisat2_index"],
    threads: config["threads"]["hisat2"]
    conda: "../../envs/align.yaml"
    log: "logs/hisat2_build.log"
    shell:
        "hisat2-build -p {threads} {input} {params.prefix} > {log} 2>&1"


rule align:
    input:
        fq1    = "results/trim/{sample}_R1.fastq.gz",
        fq2    = "results/trim/{sample}_R2.fastq.gz",
        index  = multiext(config["hisat2_index"], ".1.ht2", ".2.ht2"),
        strand = "results/strand/{sample}.strand",
    output:
        bam = "results/align/{sample}.bam",
        bai = "results/align/{sample}.bam.bai",
    params:
        index        = config["hisat2_index"],
        extra        = config["hisat2"]["extra"],
        strand_flag  = get_hisat2_strand_flag,
    threads: config["threads"]["hisat2"]
    conda: "../../envs/align.yaml"
    log: "logs/align/{sample}.log"
    shell:
        "hisat2 -x {params.index} -1 {input.fq1} -2 {input.fq2} "
        "-p {threads} {params.extra} {params.strand_flag} "
        "--summary-file {log} 2>> {log} "
        "| samtools sort -@ {threads} -o {output.bam} && "
        "samtools index {output.bam}"
