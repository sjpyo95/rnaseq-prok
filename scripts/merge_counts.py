#!/usr/bin/env python3
"""
Merge per-sample featureCounts outputs into a single gene × sample count matrix.

Usage: merge_counts.py <sample1.counts> [<sample2.counts> ...] <output.tsv>
"""
import sys
import pandas as pd

*input_files, output_file = sys.argv[1:]

dfs = []
for path in input_files:
    df = pd.read_csv(path, sep="\t", comment="#")
    # featureCounts columns: Geneid Chr Start End Strand Length <bam_path>
    bam_col  = df.columns[-1]
    sample   = bam_col.split("/")[-1].replace(".bam", "")
    dfs.append(
        df[["Geneid", bam_col]]
        .rename(columns={"Geneid": "Geneid", bam_col: sample})
        .set_index("Geneid")
    )

merged = pd.concat(dfs, axis=1)
merged.to_csv(output_file, sep="\t")
