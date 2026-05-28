# Step 4: Strandedness inference
# Samples with a known strandedness in samples.csv skip inference entirely.

rule strand_from_csv:
    """Write strandedness directly from samples.csv for known samples."""
    output:
        strand = "results/strand/{sample}.strand",
    wildcard_constraints:
        sample = "|".join(KNOWN_STRAND) if KNOWN_STRAND else "NOMATCH__",
    params:
        strand = lambda wc: samples_df.loc[wc.sample, "strandedness"],
    run:
        with open(output.strand, "w") as fh:
            fh.write(params.strand)


rule strand_subsample:
    """Subsample trimmed reads for fast alignment."""
    input:
        fq1 = "results/trim/{sample}_R1.fastq.gz",
        fq2 = "results/trim/{sample}_R2.fastq.gz",
    output:
        fq1 = temp("results/strand/{sample}_sub_R1.fastq.gz"),
        fq2 = temp("results/strand/{sample}_sub_R2.fastq.gz"),
    wildcard_constraints:
        sample = "|".join(UNKNOWN_STRAND) if UNKNOWN_STRAND else "NOMATCH__",
    params:
        n    = config["strandedness_inference"]["subsample_reads"],
        seed = 42,
    conda: "../../envs/align.yaml"
    log: "logs/strand/{sample}_subsample.log"
    shell:
        "seqtk sample -s {params.seed} {input.fq1} {params.n} | gzip > {output.fq1} 2> {log}; "
        "seqtk sample -s {params.seed} {input.fq2} {params.n} | gzip > {output.fq2} 2>> {log}"


rule strand_align:
    """Align subsampled reads and index BAM for RSeQC."""
    input:
        fq1   = "results/strand/{sample}_sub_R1.fastq.gz",
        fq2   = "results/strand/{sample}_sub_R2.fastq.gz",
        index = multiext(config["hisat2_index"], ".1.ht2", ".2.ht2"),
    output:
        bam = temp("results/strand/{sample}_sub.bam"),
        bai = temp("results/strand/{sample}_sub.bam.bai"),
    wildcard_constraints:
        sample = "|".join(UNKNOWN_STRAND) if UNKNOWN_STRAND else "NOMATCH__",
    params:
        index = config["hisat2_index"],
        extra = config["strandedness_inference"]["hisat2_extra"],
    threads: config["threads"]["hisat2"]
    conda: "../../envs/align.yaml"
    log: "logs/strand/{sample}_align.log"
    shell:
        "hisat2 -x {params.index} -1 {input.fq1} -2 {input.fq2} "
        "-p {threads} {params.extra} 2> {log} "
        "| samtools sort -@ {threads} -o {output.bam} && "
        "samtools index {output.bam}"


rule strand_infer:
    """Run RSeQC infer_experiment.py and call strand from fraction thresholds."""
    input:
        bam   = "results/strand/{sample}_sub.bam",
        bai   = "results/strand/{sample}_sub.bam.bai",
        bed12 = config["annotation_bed12"] if config["annotation_bed12"] else [],
    output:
        strand = "results/strand/{sample}.strand",
    wildcard_constraints:
        sample = "|".join(UNKNOWN_STRAND) if UNKNOWN_STRAND else "NOMATCH__",
    params:
        min_fraction = config["strandedness_inference"]["min_fraction"],
    conda: "../../envs/qc.yaml"
    log: "logs/strand/{sample}_rseqc.log"
    shell:
        "infer_experiment.py -r {input.bed12} -i {input.bam} > {log} 2>&1 && "
        "python3 scripts/parse_strand.py {log} {params.min_fraction} > {output.strand}"
