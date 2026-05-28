# CLAUDE.md

## Project

RNA-seq pipeline for prokaryotic genome analysis, built with Snakemake.
Paired-end sequencing data.

## Pipeline steps

| Step | Tool | Key options |
|------|------|-------------|
| 1. QC (raw) | FastQC | `-t {threads}` |
| 2. Trimming | fastp | `--trim_poly_g --poly_g_min_len 10 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --length_required 30` |
| 3. QC (trimmed) | FastQC + MultiQC | FastQC: `-t {threads}` / MultiQC: results/ + logs/ |
| 4. Strandedness inference | HISAT2 + RSeQC | subsample 500k reads → infer_experiment.py → min_fraction 0.7 |
| 5. Alignment | HISAT2 + samtools | `--dta --no-spliced-alignment`, `--rna-strandness FR/RF` if stranded |
| 6. Read counting | featureCounts | `-t CDS -g gene_id -p -s {0/1/2}` |
| 7. Batch correction | ComBat-seq (R) | optional; controlled by `batch_correction.enabled` in params.yaml |
| 8. DEG analysis | DESeq2 (R) | Wald test, padj < 0.05 & \|log2FC\| ≥ 1 |

## Strandedness handling

- `samples.csv` has a `strandedness` column: `FR`, `RF`, `unstranded`, or `unknown`
- Samples with `unknown` → Step 4 runs automatically
- Result is written back and used for Step 5 (`--rna-strandness`) and Step 6 (`-s`)
- featureCounts `-s` mapping: `unstranded=0`, `FR=1`, `RF=2`

## DESeq2 outputs

| File | Description |
|------|-------------|
| `deseq2_results.tsv` | Full results table |
| `deseq2_significant.tsv` | Filtered: padj < 0.05 & \|log2FC\| ≥ 1 |
| `normalized_counts.tsv` | DESeq2 size-factor normalized counts |

Plots: MA plot, PCA (raw vs normalized), heatmap (z-score of log10 normalized counts), Volcano plot

## Directory layout

| Path | Purpose |
|------|---------|
| `config/samples.csv` | Sample metadata including strandedness |
| `config/params.yaml` | All tool parameters (single source of truth) |
| `workflow/Snakefile` | Pipeline entry point |
| `workflow/rules/` | One `.smk` per step: `qc.smk`, `trim.smk`, `strand.smk`, `align.smk`, `count.smk`, `batch.smk`, `deg.smk` |
| `scripts/` | Helper Python / R scripts called by rules |
| `envs/` | Conda environment YAML files |
| `resources/` | Reference genome + GFF/GTF/BED12 (not in git) |
| `results/` | All pipeline outputs (not in git) |
| `logs/` | Per-rule logs (not in git) |

## Key conventions

- Log files: `logs/{rule}/{sample}.log`
- Outputs: `results/{step}/{sample}.*`
- Conda envs declared per-rule: `conda: "../../envs/<tool>.yaml"`
- `params.yaml` is the single source of truth for all tool flags

## Running

```bash
snakemake -n --configfile config/params.yaml          # dry run
snakemake --cores 8 --configfile config/params.yaml --use-conda
```
