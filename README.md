# RNA-seq Pipeline for Prokaryotic Genome

Snakemake-based RNA-seq analysis pipeline for prokaryotic organisms (paired-end).  
Steps 1–8: QC → Trimming → Alignment → Read counting → DESeq2.

## Getting started

### Option A — interactive setup script

```bash
git clone https://github.com/sjpyo95/rnaseq-prok.git
cd rnaseq-prok
python3 scripts/setup_wizard.py
```

The script walks you through each step and generates
`config/{experiment}/params.yaml` and `samples.csv` automatically.

### Option B — fully manual

```bash
# Copy the template and edit directly
cp -r config/example config/my_experiment
nano config/my_experiment/params.yaml   # fill in paths
nano config/my_experiment/samples.csv   # fill in samples

snakemake -n --configfile config/my_experiment/params.yaml
snakemake --cores 8 --configfile config/my_experiment/params.yaml --use-conda
```

## Pipeline overview

```
Step 1. QC — Raw reads         (FastQC)
    │
Step 2. Trimming               (fastp)
    │
Step 3. QC — Trimmed reads     (FastQC + MultiQC)
    │
Step 4. Strandedness inference (HISAT2 + RSeQC)  ← unknown samples only
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

## Configuration

Each experiment lives in its own subdirectory under `config/`:

```
config/
  example/              ← template (copy this to start)
    params.yaml         (fill in paths, set outdir: results/my_experiment)
    samples.csv         (one row per sample)
  my_experiment/        ← your config (gitignored, stays local)
    params.yaml
    samples.csv
results/
  my_experiment/        ← all outputs land here
```

### samples.csv columns

| Column | Values | Description |
|--------|--------|-------------|
| `sample_id` | string | Unique identifier |
| `condition` | string | Experimental group (e.g. `ctrl`, `treatment`) |
| `replicate` | int | Replicate number |
| `strandedness` | `FR` / `RF` / `unstranded` / `unknown` | `unknown` → auto-inferred (Step 4) |
| `fq1` | path | R1 FASTQ (absolute path) |
| `fq2` | path | R2 FASTQ (absolute path) |

### params.yaml — required fields to fill in

| Key | Description |
|-----|-------------|
| `samples` | Path to your samples.csv |
| `outdir` | Output directory (e.g. `results/my_experiment`) |
| `genome` | Reference genome FASTA |
| `annotation_gtf` | GTF annotation for featureCounts |
| `hisat2_index` | HISAT2 index prefix (build with `snakemake hisat2_build` if missing) |
| `deg.reference_condition` | Control condition name (must match samples.csv) |

## Key outputs

All files are written under `results/{experiment}/`:

| File | Description |
|------|-------------|
| `qc/multiqc_report.html` | Aggregated QC report (all steps) |
| `counts/raw_counts.tsv` | Gene × sample count matrix |
| `deg/deseq2_results.tsv` | Full DESeq2 results table |
| `deg/deseq2_significant.tsv` | Filtered: padj < 0.05 & \|log2FC\| ≥ 1 |
| `deg/normalized_counts.tsv` | Size-factor normalized counts |
| `deg/plot_MA.pdf` | MA plot |
| `deg/plot_PCA.pdf` | PCA (condition colors) |
| `deg/plot_heatmap_top25.pdf` | Top 25 DEGs per direction (z-score) |
| `deg/plot_heatmap_all.pdf` | All DEGs with Up/Down annotation bar |
| `deg/plot_volcano.pdf` | Volcano plot |

## Requirements

- Snakemake ≥ 7.0
- Conda / Mamba (per-rule environments via `--use-conda`)
- Python ≥ 3.9 (for `setup_wizard.py`; stdlib only, no extra packages)
