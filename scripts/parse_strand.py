#!/usr/bin/env python3
"""
Parse RSeQC infer_experiment.py output and call strandedness.

Usage: parse_strand.py <rseqc_log> <min_fraction>
Prints: FR | RF | unstranded
"""
import sys
import re

log_file     = sys.argv[1]
min_fraction = float(sys.argv[2])

txt = open(log_file).read()

# Paired-end output lines:
# Fraction of reads explained by "1++,1--,2+-,2-+": 0.05  → RF
# Fraction of reads explained by "1+-,1-+,2++,2--": 0.93  → FR
rf_match = re.search(r'"1\+\+,1--,2\+-,2-\+":\s*([\d.]+)', txt)
fr_match = re.search(r'"1\+-,1-\+,2\+\+,2--":\s*([\d.]+)', txt)

rf = float(rf_match.group(1)) if rf_match else 0.0
fr = float(fr_match.group(1)) if fr_match else 0.0

if fr >= min_fraction:
    print("FR")
elif rf >= min_fraction:
    print("RF")
else:
    print("unstranded")
