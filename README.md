# RNA-seq Pipeline for Prokaryotic Genome

Snakemake-based RNA-seq analysis pipeline for prokaryotic organisms (paired-end).

## Pipeline overview

```
Step 1. QC — Raw reads         (FastQC)
    │
Step 2. Trimming               (fastp)
    │
Step 3. QC — Trimmed reads     (FastQC + MultiQC)
    │
Step 4. Strandedness inference (HISAT2 + RSeQC)  ← samples with strandedness: unknown only
    │
Step 5. Alignment              (HISAT2 + samtools sort)
    │
Step 6. Read counting          (featureCounts)
    │
Step 7. Batch correction       (ComBat-seq)        ← optional
    │
Step 8. DEG analysis           (DESeq2)
         └── MA plot, PCA, Heatmap (top25 / all DEGs), Volcano plot
```

## Directory structure

```
RNA-pipeline/
├── config/
│   ├── samples.csv      # Sample sheet (sample_id, condition, replicate, strandedness, fq1, fq2)
│   └── params.yaml      # All tool parameters (single source of truth)
├── workflow/
│   ├── Snakefile        # Pipeline entry point
│   └── rules/           # One .smk file per step
├── scripts/             # Helper Python / R scripts
├── envs/                # Per-rule conda environment YAML files
├── resources/           # Reference genome + annotation (not tracked by git)
├── results/             # Pipeline outputs (not tracked by git)
├── logs/                # Execution logs (not tracked by git)
└── docs/                # Analysis notes
```

## Quick start

```bash
# 1. Edit config/samples.csv and config/params.yaml for your dataset

# 2. Dry run
snakemake -n --configfile config/params.yaml

# 3. Run (conda environments are created automatically per rule)
snakemake --cores 8 --configfile config/params.yaml --use-conda
```

## Sample sheet format (`config/samples.csv`)

| Column | Values | Description |
|--------|--------|-------------|
| `sample_id` | string | Unique sample identifier |
| `condition` | string | Experimental condition (e.g. `ctrl`, `treatment`) |
| `replicate` | int | Replicate number |
| `strandedness` | `FR` / `RF` / `unstranded` / `unknown` | `unknown` triggers automatic inference (Step 4) |
| `fq1` | path | R1 FASTQ path |
| `fq2` | path | R2 FASTQ path |

## Key outputs

| File | Description |
|------|-------------|
| `results/qc/multiqc_report.html` | Aggregated QC report |
| `results/counts/raw_counts.tsv` | Gene × sample count matrix |
| `results/deg/deseq2_results.tsv` | Full DESeq2 results |
| `results/deg/deseq2_significant.tsv` | Filtered: padj < 0.05 & \|log2FC\| ≥ 1 |
| `results/deg/normalized_counts.tsv` | DESeq2 size-factor normalized counts |
| `results/deg/plot_MA.pdf` | MA plot |
| `results/deg/plot_PCA.pdf` | PCA (condition colors) |
| `results/deg/plot_heatmap_top25.pdf` | Top 25 DEGs per direction (z-score heatmap) |
| `results/deg/plot_heatmap_all.pdf` | All DEGs with Up/Down annotation bar |
| `results/deg/plot_volcano.pdf` | Volcano plot |

## Requirements

- Snakemake >= 7.0
- Conda / Mamba (conda environments are managed per rule via `--use-conda`)
