STRAND_TO_FC = {"unstranded": "0", "FR": "1", "RF": "2"}


def get_fc_strand(wildcards):
    strand = open(f"results/strand/{wildcards.sample}.strand").read().strip()
    return STRAND_TO_FC.get(strand, "0")


rule featurecounts:
    input:
        bam    = "results/align/{sample}.bam",
        bai    = "results/align/{sample}.bam.bai",
        strand = "results/strand/{sample}.strand",
        gtf    = config["annotation_gtf"],
    output:
        counts  = "results/counts/{sample}.counts",
        summary = "results/counts/{sample}.counts.summary",
    params:
        feature   = config["featurecounts"]["feature_type"],
        attribute = config["featurecounts"]["attribute"],
        strand    = get_fc_strand,
        extra     = config["featurecounts"]["extra"],
        paired    = "-p" if config["featurecounts"]["paired"] else "",
        overlap   = "-O" if config["featurecounts"]["multi_overlap"] else "",
    threads: config["threads"]["featurecounts"]
    conda: "../../envs/count.yaml"
    log: "logs/count/{sample}.log"
    shell:
        "featureCounts "
        "-T {threads} "
        "-a {input.gtf} "
        "-t {params.feature} "
        "-g {params.attribute} "
        "-s {params.strand} "
        "{params.paired} {params.overlap} "
        "-o {output.counts} "
        "{params.extra} "
        "{input.bam} "
        "> {log} 2>&1"


rule merge_counts:
    input:
        expand("results/counts/{sample}.counts", sample=SAMPLES),
    output:
        "results/counts/raw_counts.tsv",
    conda: "../../envs/count.yaml"
    log: "logs/count/merge.log"
    shell:
        "python3 scripts/merge_counts.py {input} {output} > {log} 2>&1"
