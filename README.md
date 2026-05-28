# RNA-seq Pipeline for Prokaryotic Genome

Snakemake-based RNA-seq analysis pipeline for prokaryotic organisms.

## Pipeline overview

```
Raw reads (FASTQ)
    │
    ▼
QC (FastQC / MultiQC)
    │
    ▼
Trimming (fastp)
    │
    ▼
Alignment (HISAT2)
    │
    ▼
Read counting (featureCounts)
    │
    ▼
DEG analysis (DESeq2 / edgeR)
```

## Directory structure

```
RNA-pipeline/
├── config/          # Sample sheet and pipeline parameters
├── workflow/
│   ├── rules/       # Snakemake rules per step
│   └── Snakefile    # Main workflow entry point
├── scripts/         # Helper Python/R scripts
├── envs/            # Conda environment YAML files
├── resources/       # Reference genome and annotation (not tracked)
├── results/         # Pipeline outputs (not tracked)
├── logs/            # Execution logs (not tracked)
└── docs/            # Analysis notes and protocol
```

## Quick start

```bash
# 1. Create conda environment
conda env create -f envs/rna_pipeline.yaml
conda activate rna_pipeline

# 2. Edit config/samples.csv and config/params.yaml

# 3. Dry run
snakemake -n --configfile config/params.yaml

# 4. Run
snakemake --cores 8 --configfile config/params.yaml --use-conda
```

## Requirements

- Snakemake >= 7.0
- Conda / Mamba
