# CLAUDE.md

## Project

RNA-seq pipeline for prokaryotic genome analysis, built with Snakemake.

## Pipeline steps

1. QC — FastQC + MultiQC
2. Trimming — fastp
3. Alignment — HISAT2
4. Read counting — featureCounts
5. DEG analysis — DESeq2 / edgeR (R)

## Directory layout

| Path | Purpose |
|------|---------|
| `config/` | `samples.csv` (sample metadata), `params.yaml` (tool parameters) |
| `workflow/Snakefile` | Pipeline entry point |
| `workflow/rules/` | One `.smk` file per step (qc, trim, align, count, deg) |
| `scripts/` | Helper Python / R scripts called by rules |
| `envs/` | Conda environment YAML files |
| `resources/` | Reference genome + GFF annotation (not in git) |
| `results/` | All pipeline outputs (not in git) |
| `logs/` | Per-rule logs (not in git) |

## Key conventions

- Each Snakemake rule lives in its own file under `workflow/rules/`
- Log files go to `logs/{rule}/{sample}.log`
- Intermediate and final outputs go to `results/{step}/{sample}.*`
- Conda envs are declared per-rule with `conda: "../../envs/<tool>.yaml"`
- `config/params.yaml` is the single source of truth for all tool flags

## Running the pipeline

```bash
# Dry run
snakemake -n --configfile config/params.yaml

# Full run
snakemake --cores 8 --configfile config/params.yaml --use-conda
```
