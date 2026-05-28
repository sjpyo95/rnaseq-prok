rule trim:
    input:
        fq1 = get_fq1,
        fq2 = get_fq2,
    output:
        fq1  = "results/trim/{sample}_R1.fastq.gz",
        fq2  = "results/trim/{sample}_R2.fastq.gz",
        html = "results/trim/{sample}_fastp.html",
        json = "results/trim/{sample}_fastp.json",
    params:
        poly_g_min_len         = config["fastp"]["poly_g_min_len"],
        cut_right_window_size  = config["fastp"]["cut_right_window_size"],
        cut_right_mean_quality = config["fastp"]["cut_right_mean_quality"],
        length_required        = config["fastp"]["length_required"],
        extra                  = config["fastp"]["extra"],
    threads: config["threads"]["fastp"]
    conda: "../../envs/fastp.yaml"
    log: "logs/trim/{sample}.log"
    shell:
        "fastp "
        "-i {input.fq1} -I {input.fq2} "
        "-o {output.fq1} -O {output.fq2} "
        "-h {output.html} -j {output.json} "
        "--trim_poly_g --poly_g_min_len {params.poly_g_min_len} "
        "--cut_right "
        "--cut_right_window_size {params.cut_right_window_size} "
        "--cut_right_mean_quality {params.cut_right_mean_quality} "
        "--length_required {params.length_required} "
        "--thread {threads} "
        "{params.extra} "
        "> {log} 2>&1"
