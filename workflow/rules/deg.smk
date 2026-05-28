rule deseq2:
    input:
        counts  = COUNTS_FOR_DEG,
        samples = config["samples"],
    output:
        results     = "results/deg/deseq2_results.tsv",
        significant = "results/deg/deseq2_significant.tsv",
        norm_counts = "results/deg/normalized_counts.tsv",
        plot_ma      = "results/deg/plot_MA.pdf",
        plot_pca     = "results/deg/plot_PCA.pdf",
        plot_heatmap = "results/deg/plot_heatmap.pdf",
        plot_volcano = "results/deg/plot_volcano.pdf",
    params:
        ref_condition   = config["deg"]["reference_condition"],
        padj_cutoff     = config["deg"]["padj_cutoff"],
        lfc_cutoff      = config["deg"]["lfc_cutoff"],
        use_batch       = config["batch_correction"]["enabled"],
        batch_col       = config["batch_correction"]["batch_column"],
    conda: "../../envs/r.yaml"
    log: "logs/deg/deseq2.log"
    script:
        "../../scripts/run_deseq2.R"
