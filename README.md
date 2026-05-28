# RNA-seq Pipeline for Prokaryotic Genome

Snakemake-based RNA-seq analysis pipeline for prokaryotic organisms.

## Pipeline overview

```
Step 1. QC — Raw reads       (FastQC)
    │
Step 2. Trimming              (fastp)
    │
Step 3. QC — Trimmed reads   (FastQC + MultiQC)
    │
Step 4. Strandedness inference (HISAT2 + RSeQC) ← samples with strandedness: unknown only
    │
Step 5. Alignment             (HISAT2 + samtools sort)
    │
Step 6. Read counting         (featureCounts)
    │
Step 7. Batch correction      (ComBat-seq) ← optional
    │
Step 8. DEG analysis          (DESeq2)
         └── Outputs: MA plot, PCA, Heatmap, Volcano plot
```

## Directory structure

```
RNA-pipeline/
├── config/
│   ├── samples.csv      # Sample sheet (sample_id, condition, replicate, strandedness, fq1, fq2)
│   └── params.yaml      # All tool parameters
├── workflow/
│   ├── rules/           # One .smk file per step
│   └── Snakefile        # Pipeline entry point
├── scripts/             # Helper Python / R scripts
├── envs/                # Conda environment YAML files
├── resources/           # Reference genome + annotation (not tracked by git)
├── results/             # Pipeline outputs (not tracked by git)
├── logs/                # Execution logs (not tracked by git)
└── docs/                # Analysis notes and protocol
```

## Quick start

```bash
# 1. Create conda environment
conda env create -f envs/rna_pipeline.yaml
conda activate rna_pipeline

# 2. Fill in config/samples.csv and verify config/params.yaml

# 3. Dry run
snakemake -n --configfile config/params.yaml

# 4. Run
snakemake --cores 8 --configfile config/params.yaml --use-conda
```

## Requirements

- Snakemake >= 7.0
- Conda / Mamba
