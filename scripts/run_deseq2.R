library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(readr)
library(dplyr)
library(tibble)

log_file <- snakemake@log[[1]]
con <- file(log_file, open = "wt")
sink(con, type = "message")

# --- Colors (matching reference pipeline) ---
up_col   <- "#b2182b"
down_col <- "#2166ac"
not_col  <- "gray90"
heat_pal <- colorRampPalette(c("#2166ac", "gray90", "#b2182b"))(201)

# --- Load inputs ---
counts_mat <- read_tsv(snakemake@input[["counts"]]) |>
    column_to_rownames("Geneid") |>
    as.matrix() |>
    round()

meta <- read_csv(snakemake@input[["samples"]]) |>
    column_to_rownames("sample_id")

ref_condition <- snakemake@params[["ref_condition"]]
meta$condition <- factor(meta$condition,
                         levels = c(ref_condition,
                                    setdiff(unique(meta$condition), ref_condition)))

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

target <- setdiff(levels(meta$condition), ref_condition)[1]
padj_cutoff <- snakemake@params[["padj_cutoff"]]
lfc_cutoff  <- snakemake@params[["lfc_cutoff"]]

res <- results(dds,
               contrast = c("condition", target, ref_condition),
               alpha    = padj_cutoff)

res_df <- as.data.frame(res) |>
    rownames_to_column("gene_id") |>
    arrange(padj)

res_df$neg_log10_padj <- -log10(pmax(res_df$padj, .Machine$double.xmin))

res_df$status <- factor(
    ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange >=  lfc_cutoff, "Up-regulated",
    ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange <= -lfc_cutoff, "Down-regulated",
           "Not-significant")),
    levels = c("Down-regulated", "Not-significant", "Up-regulated")
)

sig_df <- filter(res_df,
                 status %in% c("Up-regulated", "Down-regulated"))

norm_counts <- counts(dds, normalized = TRUE) |>
    as.data.frame() |>
    rownames_to_column("gene_id")

write_tsv(res_df |> select(-status, -neg_log10_padj), snakemake@output[["results"]])
write_tsv(sig_df |> select(-status, -neg_log10_padj), snakemake@output[["significant"]])
write_tsv(norm_counts,               snakemake@output[["norm_counts"]])

# --- Prepare count matrix (control first, target second) ---
control_samples <- rownames(meta)[meta$condition == ref_condition]
target_samples  <- rownames(meta)[meta$condition == target]
sample_order    <- c(control_samples, target_samples)

norm_mat <- norm_counts |> column_to_rownames("gene_id") |> as.matrix()
norm_mat <- norm_mat[, sample_order, drop = FALSE]

ann_col <- data.frame(
    condition = factor(
        c(rep(ref_condition, length(control_samples)),
          rep(target,         length(target_samples))),
        levels = c(ref_condition, target)
    ),
    row.names = sample_order
)
ann_col_colors <- list(
    condition = setNames(c(down_col, up_col), c(ref_condition, target))
)

# --- Sort DEGs: down (by |lfc|) then up (by |lfc|) ---
downs_all <- sig_df[sig_df$status == "Down-regulated", ]
ups_all   <- sig_df[sig_df$status == "Up-regulated",   ]
downs_all <- downs_all[order(downs_all$log2FoldChange),   , drop = FALSE]
ups_all   <- ups_all[order(-ups_all$log2FoldChange),      , drop = FALSE]

# --- MA plot ---
pdf(snakemake@output[["plot_ma"]])
plotMA(res, alpha = padj_cutoff)
dev.off()

# --- PCA ---
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))
cond_colors <- setNames(c(down_col, up_col), c(ref_condition, target))
pdf(snakemake@output[["plot_pca"]])
print(
    ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
        geom_point(size = 3) +
        geom_text_repel(size = 3, show.legend = FALSE) +
        scale_color_manual(values = cond_colors) +
        xlab(paste0("PC1: ", pct_var[1], "% variance")) +
        ylab(paste0("PC2: ", pct_var[2], "% variance")) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(sprintf("PCA - %s vs %s", target, ref_condition))
)
dev.off()

# --- Heatmap helper ---
make_heatmap <- function(gene_ids, show_rownames, ann_row = NULL, ann_row_colors = NULL,
                         title, outfile) {
    gene_ids <- gene_ids[gene_ids %in% rownames(norm_mat)]
    if (length(gene_ids) < 2) {
        pdf(outfile); plot.new(); text(0.5, 0.5, "Not enough DEGs"); dev.off()
        return(invisible(NULL))
    }
    mat <- log2(norm_mat[gene_ids, , drop = FALSE] + 1)
    ann_colors <- ann_col_colors
    if (!is.null(ann_row_colors)) ann_colors <- c(ann_colors, ann_row_colors)
    n_genes   <- length(gene_ids)
    fig_height <- max(6, 0.18 * n_genes + 3)
    # width: enough for heatmap body + row labels + legend + title headroom
    fig_width  <- max(10, ncol(mat) * 0.7 + 6)
    pdf(outfile, width = fig_width, height = fig_height)
    pheatmap(
        mat,
        color             = heat_pal,
        scale             = "row",
        cluster_rows      = FALSE,
        cluster_cols      = FALSE,
        annotation_col    = ann_col,
        annotation_row    = ann_row,
        annotation_colors = ann_colors,
        show_rownames     = show_rownames,
        show_colnames     = TRUE,
        cellwidth         = 14,
        main              = " "    # reserve space; real title drawn below
    )
    # Overlay title in root viewport (always page-centred).
    # pheatmap reserves 1.5 cm for main=" "; match that height in npc.
    title_npc <- as.numeric(grid::convertHeight(grid::unit(1.5, "cm"), "npc"))
    grid::upViewport(0)
    grid::pushViewport(grid::viewport(y = 1, height = title_npc, just = "top"))
    grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
    grid::grid.text(title, x = 0.5, y = 0.5,
                    gp = grid::gpar(fontface = "bold", fontsize = 14))
    grid::popViewport()
    dev.off()
}

# Top-25 heatmap (top 25 per direction by |log2FC|)
top25_ids <- c(head(downs_all$gene_id, 25), head(ups_all$gene_id, 25))
n_down25  <- min(25, nrow(downs_all))
n_up25    <- min(25, nrow(ups_all))
top25_title <- sprintf("%s vs %s Top 25 DEGs", target, ref_condition)
make_heatmap(top25_ids, show_rownames = TRUE, title = top25_title,
             outfile = snakemake@output[["plot_heatmap_top25"]])

# All-DEGs heatmap (with left row annotation)
all_ids  <- c(downs_all$gene_id, ups_all$gene_id)
n_down   <- nrow(downs_all)
n_up     <- nrow(ups_all)
if (length(all_ids) >= 2) {
    ann_row <- data.frame(
        status = factor(
            c(rep("Down-regulated", n_down), rep("Up-regulated", n_up)),
            levels = c("Down-regulated", "Up-regulated")
        ),
        row.names = all_ids
    )
    ann_row_colors <- list(
        status = c("Down-regulated" = down_col, "Up-regulated" = up_col)
    )
    all_title <- sprintf("%s vs %s Top all DEGs (Total = %d)", target, ref_condition, n_down + n_up)
    make_heatmap(all_ids, show_rownames = FALSE,
                 ann_row = ann_row, ann_row_colors = ann_row_colors,
                 title = all_title, outfile = snakemake@output[["plot_heatmap_all"]])
} else {
    pdf(snakemake@output[["plot_heatmap_all"]])
    plot.new(); text(0.5, 0.5, "Not enough DEGs for heatmap")
    dev.off()
}

# --- Volcano plot ---
res_plot <- res_df
upLab   <- head(ups_all[order(-ups_all$log2FoldChange), ], 15)
downLab <- head(downs_all[order(downs_all$log2FoldChange), ], 15)
label_df <- rbind(downLab, upLab)
max_x    <- max(abs(res_plot$log2FoldChange), na.rm = TRUE)
finite_y <- res_plot$neg_log10_padj[is.finite(res_plot$neg_log10_padj)]
y_cap <- if (length(finite_y) > 0) max(finite_y, na.rm = TRUE) else 10

pdf(snakemake@output[["plot_volcano"]], width = 6, height = 6)
print(
    ggplot(res_plot, aes(x = log2FoldChange, y = neg_log10_padj, color = status)) +
        geom_point(alpha = 0.8, size = 1.5) +
        scale_color_manual(values = c(
            "Down-regulated"  = down_col,
            "Not-significant" = not_col,
            "Up-regulated"    = up_col
        )) +
        geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "black") +
        geom_hline(yintercept = -log10(padj_cutoff),         linetype = "dashed", color = "black") +
        scale_x_continuous(limits = c(-max_x, max_x)) +
        scale_y_continuous(limits = c(0, y_cap)) +
        geom_text_repel(
            data          = label_df,
            aes(label     = gene_id),
            size          = 2.5,
            max.overlaps  = Inf,
            box.padding   = 0.3,
            point.padding = 0.2,
            min.segment.length = 0,
            show.legend   = FALSE
        ) +
        labs(
            title = sprintf("Volcano Plot of %s vs %s", target, ref_condition),
            x     = "log2(Fold Change)",
            y     = "-log10(Adjusted P-value)"
        ) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank())
)
dev.off()

sink(type = "message")
close(con)
