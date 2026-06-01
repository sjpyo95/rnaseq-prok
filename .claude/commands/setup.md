Guide the user through first-time setup of this RNA-seq pipeline.
Run this whenever: no experiment config exists yet, or the user wants to add a new experiment.

Follow these steps in order. Be friendly and concise. Use AskUserQuestion for user decisions.
Run Bash commands to validate paths and generate files. Show progress clearly.

---

## Step 1: Welcome

Greet the user. Explain what will happen:
- We'll create a config for their experiment (takes ~5 minutes)
- They need: raw FASTQ files, reference genome (.fa), GTF annotation, HISAT2 index
- At the end they'll have a ready-to-run command

Check if any experiments already exist:
```bash
find config -name "params.yaml" -mindepth 2 2>/dev/null
```
If experiments exist, list them and ask whether to add a new one or work with an existing one.

---

## Step 2: Experiment name

Ask the user for a short experiment name.
- Examples: `heat_stress`, `sodium_vs_ctrl`, `anaerobic`
- Rules: no spaces, lowercase, underscores OK
- This becomes the output directory: `results/<name>/`

---

## Step 3: Scan for raw FASTQ files

Ask for the raw data directory path.

Run:
```bash
python3 scripts/scan_fastq.py <data_dir>
```

Show the discovered samples as a table. If no files found:
- Check if the path exists (`ls <path>`)
- Try the parent directory
- Ask if files have a different naming convention

If R2 files are missing for any sample, flag this clearly and ask the user to verify.

---

## Step 4: Assign conditions and replicates

Show the draft sample table (sample_id, fq1, fq2).

Ask:
1. How many experimental conditions? (e.g. 2 — control and treatment)
2. Which samples belong to which condition? Let the user match by sample name.
3. What is the control/reference condition name? (exact name used for DESeq2)
4. Replicate numbers for each sample (1, 2, 3…)

Build the final samples.csv table and show it for confirmation before writing.

---

## Step 5: Strandedness

Ask: "Do you know the library strandedness?"

Options to present:
- **RF** — most common (Illumina TruSeq Stranded mRNA, dUTP method)
- **FR** — less common
- **unstranded** — non-directional library prep
- **unknown** — pipeline will infer automatically (requires BED12 annotation file)

Guidance if unsure:
- TruSeq Stranded → RF
- TruSeq (non-stranded) → unstranded
- If unknown: ask if they have a BED12 file for the genome. If yes → unknown. If no → ask their kit.

Set all samples to the same strandedness value (or unknown to auto-detect).

---

## Step 6: Reference genome files

Ask for each path, then validate with `ls -lh <path>`:

1. **Reference genome FASTA** — check file exists and is non-empty
2. **GTF annotation** — check file exists
3. **HISAT2 index prefix** — check that `<prefix>.1.ht2` exists:
   ```bash
   ls <prefix>.1.ht2
   ```

If the HISAT2 index doesn't exist:
- Explain it needs to be built once from the genome (~10–30 min)
- Ask: "Should I build the index now, or will you provide a pre-built one?"
- If build now: note the command they'll need to run:
  ```bash
  snakemake hisat2_build --cores 8 --configfile config/<exp>/params.yaml --use-conda
  ```
  (They should run this before the full pipeline)

---

## Step 7: Generate config files

Create `config/<experiment>/` with two files.

**samples.csv** — use the table built in Steps 4–5:
```
sample_id,condition,replicate,strandedness,fq1,fq2
```

**params.yaml** — copy `config/awrp_sodium/params.yaml` as a template, then update:
- `samples: config/<experiment>/samples.csv`
- `outdir: results/<experiment>`
- `genome:`, `annotation_gtf:`, `hisat2_index:` from Step 6
- `deg.reference_condition:` from Step 4
- Leave all other tool parameters at defaults unless the user asks to change them

Show both files after writing.

---

## Step 8: Dry-run validation

Run:
```bash
snakemake -n --configfile config/<experiment>/params.yaml 2>&1 | head -50
```

If it succeeds: show the job count summary (N jobs, rules breakdown).

If it fails, diagnose:
- `KeyError` in config → missing key in params.yaml → fix and re-run
- Missing input files → paths in samples.csv or params.yaml are wrong → check and fix
- Strandedness wildcard error → wildcard_constraints mismatch → check samples.csv

---

## Step 9: Ready

Show a final summary:
```
Experiment : <name>
Samples    : N (X conditions)
Output dir : results/<experiment>/
Config     : config/<experiment>/params.yaml
```

Show the run command:
```bash
snakemake --cores 8 --configfile config/<experiment>/params.yaml --use-conda
```

Ask: "Would you like to start the pipeline now, or run it later?"
- If now: run it
- If later: remind them to use the command above
