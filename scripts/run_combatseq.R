library(sva)
library(DESeq2)
library(readr)
library(dplyr)

log_file <- snakemake@log[[1]]
con <- file(log_file, open = "wt")
sink(con, type = "message")

counts_mat <- read_tsv(snakemake@input[["counts"]]) |>
    column_to_rownames("Geneid") |>
    as.matrix() |>
    round()

meta <- read_csv(snakemake@input[["samples"]])
batch_col <- snakemake@params[["batch_col"]]

stopifnot(batch_col %in% colnames(meta))
meta <- meta[match(colnames(counts_mat), meta[["sample_id"]]), ]
stopifnot(!anyNA(meta[["sample_id"]]))
batch     <- meta[[batch_col]]
condition <- meta[["condition"]]

corrected <- ComBat_seq(counts_mat, batch = batch, group = condition)

# VST for QC plots (written to log for inspection, not saved as output)
dds_vst <- DESeqDataSetFromMatrix(corrected, colData = meta, design = ~condition)
vsd     <- vst(dds_vst, blind = TRUE)
message("VST complete — use plotPCA(vsd) for batch QC inspection")

write_tsv(as.data.frame(corrected) |> tibble::rownames_to_column("Geneid"),
          snakemake@output[["corrected"]])

sink(type = "message")
close(con)
