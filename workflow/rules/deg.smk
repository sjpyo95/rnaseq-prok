rule deseq2:
    input:
        counts  = COUNTS_FOR_DEG,
        samples = config["samples"],
    output:
        results            = OUTDIR + "/deg/deseq2_results.tsv",
        significant        = OUTDIR + "/deg/deseq2_significant.tsv",
        norm_counts        = OUTDIR + "/deg/normalized_counts.tsv",
        plot_ma            = OUTDIR + "/deg/plot_MA.pdf",
        plot_pca           = OUTDIR + "/deg/plot_PCA.pdf",
        plot_heatmap_top25 = OUTDIR + "/deg/plot_heatmap_top25.pdf",
        plot_heatmap_all   = OUTDIR + "/deg/plot_heatmap_all.pdf",
        plot_volcano       = OUTDIR + "/deg/plot_volcano.pdf",
    params:
        ref_condition = config["deg"]["reference_condition"],
        padj_cutoff   = config["deg"]["padj_cutoff"],
        lfc_cutoff    = config["deg"]["lfc_cutoff"],
        use_batch     = config["batch_correction"]["enabled"],
        batch_col     = config["batch_correction"]["batch_column"],
    conda: "../../envs/r.yaml"
    log: "logs/deg/deseq2.log"
    script:
        "../../scripts/run_deseq2.R"
