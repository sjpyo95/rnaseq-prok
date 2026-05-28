library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(ggplot2)
library(readr)
library(dplyr)
library(tibble)

log_file <- snakemake@log[[1]]
con <- file(log_file, open = "wt")
sink(con, type = "message")

# --- Load inputs ---
counts_mat <- read_tsv(snakemake@input[["counts"]]) |>
    column_to_rownames("Geneid") |>
    as.matrix() |>
    round()

meta <- read_csv(snakemake@input[["samples"]]) |>
    column_to_rownames("sample_id")
meta$condition <- factor(meta$condition,
                         levels = c(snakemake@params[["ref_condition"]],
                                    setdiff(unique(meta$condition),
                                            snakemake@params[["ref_condition"]])))

# --- Design formula ---
design_formula <- if (snakemake@params[["use_batch"]]) {
    as.formula(paste0("~ ", snakemake@params[["batch_col"]], " + condition"))
} else {
    ~ condition
}

# --- DESeq2 ---
dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData   = meta,
                              design    = design_formula)
dds <- DESeq(dds)

target <- setdiff(levels(meta$condition), snakemake@params[["ref_condition"]])[1]
res    <- results(dds,
                  contrast = c("condition", target, snakemake@params[["ref_condition"]]),
                  alpha    = snakemake@params[["padj_cutoff"]])

res_df  <- as.data.frame(res) |> rownames_to_column("gene_id") |> arrange(padj)
sig_df  <- filter(res_df,
                  padj < snakemake@params[["padj_cutoff"]],
                  abs(log2FoldChange) >= snakemake@params[["lfc_cutoff"]])

norm_counts <- counts(dds, normalized = TRUE) |>
    as.data.frame() |>
    rownames_to_column("gene_id")

write_tsv(res_df,     snakemake@output[["results"]])
write_tsv(sig_df,     snakemake@output[["significant"]])
write_tsv(norm_counts, snakemake@output[["norm_counts"]])

# --- Plots ---
vsd <- vst(dds, blind = FALSE)

# MA plot
pdf(snakemake@output[["plot_ma"]])
plotMA(res, alpha = snakemake@params[["padj_cutoff"]])
dev.off()

# PCA
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))
pdf(snakemake@output[["plot_pca"]])
ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5, size = 3) +
    xlab(paste0("PC1: ", pct_var[1], "% variance")) +
    ylab(paste0("PC2: ", pct_var[2], "% variance")) +
    theme_bw()
dev.off()

# Heatmap (z-score of log10 normalized counts, top 50 DEGs)
top_genes <- head(sig_df$gene_id, 50)
if (length(top_genes) > 1) {
    mat_all <- norm_counts |> column_to_rownames("gene_id") |> as.matrix()
    mat <- log10(mat_all[top_genes, ] + 1)
    mat_z <- t(scale(t(mat)))
    pdf(snakemake@output[["plot_heatmap"]])
    pheatmap(mat_z,
             annotation_col = meta["condition"],
             show_rownames   = length(top_genes) <= 30,
             fontsize_row    = 7,
             main            = "Top DEGs — z-score log10(norm counts)")
    dev.off()
} else {
    pdf(snakemake@output[["plot_heatmap"]])
    plot.new(); text(0.5, 0.5, "Not enough significant DEGs for heatmap")
    dev.off()
}

# Volcano plot
pdf(snakemake@output[["plot_volcano"]])
EnhancedVolcano(res_df,
    lab      = res_df$gene_id,
    x        = "log2FoldChange",
    y        = "padj",
    pCutoff  = snakemake@params[["padj_cutoff"]],
    FCcutoff = snakemake@params[["lfc_cutoff"]],
    title    = paste(target, "vs", snakemake@params[["ref_condition"]]))
dev.off()

sink(type = "message")
close(con)
