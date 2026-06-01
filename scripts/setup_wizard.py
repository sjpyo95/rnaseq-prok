#!/usr/bin/env python3
"""
Interactive setup wizard for the RNA-seq pipeline.
Creates config/{experiment}/params.yaml and samples.csv.

Usage:
    python3 scripts/setup_wizard.py
"""

import sys
import re
import subprocess
from pathlib import Path

# ── ANSI colors (gracefully disabled on non-TTY) ──────────────────────────────
def _ansi(code):
    return code if sys.stdout.isatty() else ""

B   = _ansi("\033[1m");  DIM = _ansi("\033[2m");  RST = _ansi("\033[0m")
G   = _ansi("\033[32m"); C   = _ansi("\033[36m")
Y   = _ansi("\033[33m"); R   = _ansi("\033[31m")

def hdr(text):
    print(f"\n{B}{C}{'─' * 58}{RST}")
    print(f"{B}{C}  {text}{RST}")
    print(f"{B}{C}{'─' * 58}{RST}")

def ok(text):   print(f"  {G}✓{RST}  {text}")
def warn(text): print(f"  {Y}!{RST}  {text}")
def err(text):  print(f"  {R}✗{RST}  {text}")
def info(text): print(f"     {DIM}{text}{RST}")


# ── Input helpers ─────────────────────────────────────────────────────────────
def ask(prompt, default=None, validator=None):
    """Prompt the user until a valid answer is given."""
    hint = f" {DIM}[{default}]{RST}" if default is not None else ""
    while True:
        try:
            ans = input(f"\n  {B}>{RST} {prompt}{hint}\n    ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\n\nAborted.")
            sys.exit(0)
        if not ans and default is not None:
            ans = str(default)
        if not ans:
            err("Please enter a value.")
            continue
        if validator:
            msg = validator(ans)
            if msg is not True:
                err(msg)
                continue
        return ans


def ask_path(prompt, must_exist=True, is_dir=False):
    def v(s):
        p = Path(s).expanduser().resolve()
        if must_exist and not p.exists():
            return f"Not found: {p}"
        if is_dir and p.exists() and not p.is_dir():
            return f"Not a directory: {p}"
        return True
    raw = ask(prompt, validator=v)
    return Path(raw).expanduser().resolve()


def ask_yn(prompt, default="y"):
    opts = f"{B}Y{RST}/n" if default == "y" else f"y/{B}N{RST}"
    ans = ask(f"{prompt} ({opts})", default=default)
    return ans.lower().startswith("y")


# ── FASTQ scanning ────────────────────────────────────────────────────────────
def find_r1_files(root):
    patterns = [
        "*_R1.fastq.gz", "*_R1.fq.gz",
        "*_R1_001.fastq.gz", "*_R1_001.fq.gz",
        "*_1.fastq.gz",  "*_1.fq.gz",
    ]
    found = []
    for p in patterns:
        found.extend(root.rglob(p))
    return sorted(set(found))


def infer_r2(r1):
    name = r1.name
    for src, dst in [("_R1_001.", "_R2_001."), ("_R1.", "_R2."), ("_1.", "_2.")]:
        if src in name:
            r2 = r1.parent / name.replace(src, dst, 1)
            if r2.exists():
                return r2
    return None


def infer_sample_id(r1):
    name = r1.name
    name = re.sub(r"_R1(_001)?\.(fastq|fq)\.gz$", "", name)
    name = re.sub(r"_1\.(fastq|fq)\.gz$", "", name)
    return name


# ── Config file writers (no PyYAML dependency) ────────────────────────────────
def write_params(path, samples_csv, outdir, genome, gtf, hisat2_index,
                 threads, ref_cond, padj, lfc):
    path.write_text(f"""\
samples: {samples_csv}
outdir:  {outdir}

genome:           {genome}
annotation_gtf:   {gtf}
annotation_gff:   ""
annotation_bed12: ""

hisat2_index: {hisat2_index}

threads:
  fastqc: {threads}
  fastp: {threads}
  hisat2: {threads}
  samtools: {threads}
  featurecounts: {threads}

fastqc:
  extra: ""

fastp:
  trim_poly_g: true
  poly_g_min_len: 10
  cut_right: true
  cut_right_window_size: 4
  cut_right_mean_quality: 20
  length_required: 30
  extra: ""

strandedness_inference:
  subsample_reads: 500000
  min_fraction: 0.7
  hisat2_extra: "--dta --no-spliced-alignment"

hisat2:
  extra: "--dta --no-spliced-alignment"

featurecounts:
  feature_type: "CDS"
  attribute: "gene_id"
  paired: true
  multi_overlap: false
  extra: ""

batch_correction:
  enabled: false
  batch_column: "batch"

deg:
  reference_condition: "{ref_cond}"
  padj_cutoff: {padj}
  lfc_cutoff: {lfc}
""")


# ── Main wizard ───────────────────────────────────────────────────────────────
def main():
    print(f"\n{B}{'═' * 58}{RST}")
    print(f"{B}  RNA-seq Pipeline — Setup Wizard{RST}")
    print(f"{B}{'═' * 58}{RST}")
    print("  Guides you through creating an experiment config.")
    print("  Press Ctrl-C at any time to abort.\n")

    # ── Step 1: Experiment name ───────────────────────────────────────────────
    hdr("Step 1 / 7  ·  Experiment name")
    info("This becomes the config directory and output prefix.")
    info("Example: heat_stress  →  config/heat_stress/  results/heat_stress/")

    def valid_name(s):
        if " " in s:
            return "No spaces — use underscores (e.g. heat_stress)"
        if not re.match(r"^[a-zA-Z0-9_\-]+$", s):
            return "Use only letters, numbers, underscores, or hyphens"
        if (Path("config") / s).exists():
            return f"config/{s}/ already exists — choose a different name"
        return True

    exp = ask("Experiment name", validator=valid_name)
    ok(f"Experiment: {exp}")

    # ── Step 2: Scan FASTQs ───────────────────────────────────────────────────
    hdr("Step 2 / 7  ·  Raw FASTQ files")
    info("Subdirectories are searched recursively.")
    info("Expected naming: sample_R1.fastq.gz / sample_R2.fastq.gz")

    while True:
        data_dir = ask_path("Raw data directory", must_exist=True, is_dir=True)
        r1_files = find_r1_files(data_dir)
        if not r1_files:
            err(f"No paired FASTQ files found in {data_dir}")
            if not ask_yn("Try a different path?"):
                sys.exit(1)
            continue
        break

    samples = []
    for r1 in r1_files:
        r2  = infer_r2(r1)
        sid = infer_sample_id(r1)
        samples.append({"id": sid, "fq1": r1, "fq2": r2})

    print(f"\n  {G}{B}Found {len(samples)} sample(s):{RST}")
    for i, s in enumerate(samples, 1):
        status = f"{G}R2 ✓{RST}" if s["fq2"] else f"{R}R2 missing{RST}"
        print(f"    {B}{i:2d}.{RST}  {s['id']}  {DIM}({status}{DIM}){RST}")

    missing = [s for s in samples if not s["fq2"]]
    if missing:
        warn(f"{len(missing)} sample(s) missing R2 — pipeline requires paired-end reads.")
        if not ask_yn("Continue anyway?", default="n"):
            sys.exit(1)

    # ── Step 3: Conditions & replicates ──────────────────────────────────────
    hdr("Step 3 / 7  ·  Conditions & replicates")
    info("Tell us which samples belong to which experimental group.")

    n_samples = len(samples)

    def valid_ncond(s):
        if not s.isdigit():
            return "Enter a number"
        if not (1 <= int(s) <= n_samples):
            return f"Must be between 1 and {n_samples}"
        return True

    n_cond = int(ask(
        f"How many conditions?  (e.g. 2 for ctrl + treatment)",
        validator=valid_ncond,
    ))

    assigned  = {}   # sample_id → {"condition": str, "replicate": int}
    conditions = []

    for c in range(1, n_cond + 1):
        cname = ask(f"Condition {c} name  (e.g. ctrl, heat_stress)")
        conditions.append(cname)

        unassigned = [s for s in samples if s["id"] not in assigned]
        if not unassigned:
            break

        if c == n_cond:
            # Last condition: remaining samples assigned automatically
            chosen = unassigned
            print(f"\n  Remaining samples → '{cname}':")
        else:
            print(f"\n  Available samples:")
            for s in unassigned:
                print(f"    {B}{samples.index(s) + 1:2d}.{RST}  {s['id']}")

            def valid_sel(s):
                try:
                    nums = [int(x.strip()) for x in s.split(",")]
                except ValueError:
                    return "Comma-separated numbers (e.g. 1,3)"
                for n in nums:
                    if not (1 <= n <= n_samples):
                        return f"Number out of range: {n}"
                    if samples[n - 1]["id"] in assigned:
                        return f"Sample {n} already assigned"
                return True

            sel_str = ask(
                f"Samples for '{cname}'  (comma-separated numbers)",
                validator=valid_sel,
            )
            nums   = [int(x.strip()) for x in sel_str.split(",")]
            chosen = [samples[n - 1] for n in nums]

        for rep, s in enumerate(chosen, 1):
            assigned[s["id"]] = {"condition": cname, "replicate": rep}
        ok(f"'{cname}': {len(chosen)} sample(s), replicates 1–{len(chosen)}")

    # Assign any stragglers to the last condition
    for s in samples:
        if s["id"] not in assigned:
            assigned[s["id"]] = {"condition": conditions[-1], "replicate": 1}

    print(f"\n  Conditions defined: {', '.join(conditions)}")
    ref_cond = ask(
        "Which condition is the control/reference?  (DESeq2 baseline)",
        validator=lambda s: True if s in conditions else f"Must be one of: {', '.join(conditions)}",
    )
    ok(f"Reference condition: {ref_cond}")

    # ── Step 4: Strandedness ──────────────────────────────────────────────────
    hdr("Step 4 / 7  ·  Library strandedness")
    info("TruSeq Stranded mRNA (dUTP)  →  RF  (most common)")
    info("TruSeq non-stranded          →  unstranded")
    info("Not sure                     →  unknown  (auto-detected; needs BED12)")

    STRAND_OPTS = {"RF", "FR", "unstranded", "unknown"}
    strand = ask(
        "Strandedness  [RF / FR / unstranded / unknown]",
        default="RF",
        validator=lambda s: True if s in STRAND_OPTS else f"Choose from: {', '.join(sorted(STRAND_OPTS))}",
    )
    ok(f"Strandedness: {strand}")

    # ── Step 5: Reference files ───────────────────────────────────────────────
    hdr("Step 5 / 7  ·  Reference genome files")

    genome = ask_path("Reference genome FASTA  (.fa or .fasta)")
    ok(f"Genome: {genome}")

    gtf = ask_path("GTF annotation file  (used by featureCounts)")
    ok(f"GTF: {gtf}")

    info("Enter the HISAT2 index prefix (not the .ht2 file itself)")
    info("Example: /data/ref/hisat2/genome  (genome.1.ht2 must exist there)")
    while True:
        idx_prefix = Path(ask("HISAT2 index prefix")).expanduser().resolve()
        if Path(str(idx_prefix) + ".1.ht2").exists():
            ok(f"HISAT2 index: {idx_prefix}")
            break
        warn(f"{idx_prefix}.1.ht2 not found.")
        info(f"Build it later:  snakemake hisat2_build --configfile config/{exp}/params.yaml --use-conda")
        if ask_yn("Use this path anyway and build the index later?", default="y"):
            break

    # ── Step 6: Run parameters ────────────────────────────────────────────────
    hdr("Step 6 / 7  ·  Run parameters")

    threads = int(ask(
        "CPU threads per rule",
        default="8",
        validator=lambda s: True if s.isdigit() and int(s) > 0 else "Enter a positive integer",
    ))
    padj = ask("DESeq2 adjusted p-value cutoff", default="0.05")
    lfc  = ask("DESeq2 |log2 fold-change| cutoff", default="1.0")

    # ── Step 7: Write config files ────────────────────────────────────────────
    hdr("Step 7 / 7  ·  Generating config files")

    config_dir  = Path("config") / exp
    config_dir.mkdir(parents=True, exist_ok=True)
    csv_path    = config_dir / "samples.csv"
    params_path = config_dir / "params.yaml"

    # samples.csv
    rows = ["sample_id,condition,replicate,strandedness,fq1,fq2"]
    for s in samples:
        meta = assigned[s["id"]]
        fq2  = str(s["fq2"]) if s["fq2"] else "MISSING"
        rows.append(f"{s['id']},{meta['condition']},{meta['replicate']},{strand},{s['fq1']},{fq2}")
    csv_path.write_text("\n".join(rows) + "\n")
    ok(f"Wrote  {csv_path}")

    # params.yaml
    write_params(
        path=params_path,
        samples_csv=str(csv_path),
        outdir=f"results/{exp}",
        genome=str(genome),
        gtf=str(gtf),
        hisat2_index=str(idx_prefix),
        threads=threads,
        ref_cond=ref_cond,
        padj=padj,
        lfc=lfc,
    )
    ok(f"Wrote  {params_path}")

    # ── Dry-run validation ────────────────────────────────────────────────────
    print()
    if ask_yn("Run a dry-run to validate the config now?", default="y"):
        try:
            result = subprocess.run(
                ["snakemake", "-n", "--configfile", str(params_path)],
                capture_output=True, text=True,
            )
            if result.returncode == 0:
                for line in result.stdout.splitlines():
                    if any(k in line.lower() for k in ("total", "job stats", "count")):
                        info(line.strip())
                ok("Dry-run passed — pipeline is ready!")
            else:
                err("Dry-run failed. Last error lines:")
                for line in (result.stderr or result.stdout).splitlines()[-15:]:
                    print(f"     {R}{line}{RST}")
                warn("Fix the errors above, then re-run manually.")
        except FileNotFoundError:
            warn("'snakemake' not found in PATH — skipping dry-run.")
            info("Run manually: snakemake -n --configfile " + str(params_path))

    # ── Done ──────────────────────────────────────────────────────────────────
    print(f"\n{B}{'═' * 58}{RST}")
    print(f"{B}{G}  Setup complete!{RST}")
    print(f"{B}{'═' * 58}{RST}\n")
    print(f"  Run command:")
    print(f"\n  {B}{C}snakemake --cores {threads} \\{RST}")
    print(f"  {B}{C}    --configfile config/{exp}/params.yaml \\{RST}")
    print(f"  {B}{C}    --use-conda{RST}\n")


if __name__ == "__main__":
    main()
