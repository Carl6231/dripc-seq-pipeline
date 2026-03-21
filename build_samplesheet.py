#!/usr/bin/env python3
"""
build_samplesheet.py
====================
Scan a FASTQ directory and auto-generate config/samples.tsv for the
DRIPc-seq pipeline.

Supports two naming conventions:
  1) IP-only:   {group}{rep}_R1.fastq.gz / {group}{rep}_R2.fastq.gz
  2) IP+input:  {group}{rep}_IP_R1.fastq.gz   + {group}{rep}_input_R1.fastq.gz

Usage:
  python scripts/build_samplesheet.py \\
      --fastq-dir /path/to/raw_fastq \\
      -o config/samples.tsv \\
      --groups Sham,3day,7day,14day,28day \\
      --reps 3
"""

import argparse
import os
import re
import sys


def find_fastqs(fastq_dir, groups, reps):
    """
    Detect FASTQ files and return a list of sample dicts.
    """
    samples = []
    files = sorted(os.listdir(fastq_dir))

    for grp in groups:
        for rep in range(1, reps + 1):
            sid = f"{grp}{rep}"

            # --- try IP+input naming first ---
            ip_r1 = f"{sid}_IP_R1.fastq.gz"
            ip_r2 = f"{sid}_IP_R2.fastq.gz"
            inp_r1 = f"{sid}_input_R1.fastq.gz"
            inp_r2 = f"{sid}_input_R2.fastq.gz"

            if ip_r1 in files and ip_r2 in files:
                has_input = inp_r1 in files and inp_r2 in files
                samples.append({
                    "sample_id": sid,
                    "group": grp,
                    "rep": str(rep),
                    "ip_fq1": os.path.join(fastq_dir, ip_r1),
                    "ip_fq2": os.path.join(fastq_dir, ip_r2),
                    "input_fq1": os.path.join(fastq_dir, inp_r1) if has_input else "NA",
                    "input_fq2": os.path.join(fastq_dir, inp_r2) if has_input else "NA",
                    "has_input": str(has_input).lower(),
                })
                continue

            # --- try simple naming ---
            simple_r1 = f"{sid}_R1.fastq.gz"
            simple_r2 = f"{sid}_R2.fastq.gz"

            if simple_r1 in files and simple_r2 in files:
                samples.append({
                    "sample_id": sid,
                    "group": grp,
                    "rep": str(rep),
                    "ip_fq1": os.path.join(fastq_dir, simple_r1),
                    "ip_fq2": os.path.join(fastq_dir, simple_r2),
                    "input_fq1": "NA",
                    "input_fq2": "NA",
                    "has_input": "false",
                })
                continue

            # --- also try with underscore between group and rep ---
            alt_r1 = f"{grp}_{rep}_R1.fastq.gz"
            alt_r2 = f"{grp}_{rep}_R2.fastq.gz"
            if alt_r1 in files and alt_r2 in files:
                # check for input
                alt_inp_r1 = f"{grp}_{rep}_input_R1.fastq.gz"
                alt_inp_r2 = f"{grp}_{rep}_input_R2.fastq.gz"
                has_input = alt_inp_r1 in files and alt_inp_r2 in files
                samples.append({
                    "sample_id": sid,
                    "group": grp,
                    "rep": str(rep),
                    "ip_fq1": os.path.join(fastq_dir, alt_r1),
                    "ip_fq2": os.path.join(fastq_dir, alt_r2),
                    "input_fq1": os.path.join(fastq_dir, alt_inp_r1) if has_input else "NA",
                    "input_fq2": os.path.join(fastq_dir, alt_inp_r2) if has_input else "NA",
                    "has_input": str(has_input).lower(),
                })
                continue

            print(f"WARNING: Could not find FASTQ files for sample {sid}", file=sys.stderr)

    return samples


def write_tsv(samples, output):
    """Write the sample sheet."""
    header = ["sample_id", "group", "rep", "ip_fq1", "ip_fq2",
              "input_fq1", "input_fq2", "has_input"]
    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
    with open(output, "w") as fh:
        fh.write("# Auto-generated sample sheet for DRIPc-seq pipeline\n")
        fh.write("\t".join(header) + "\n")
        for s in samples:
            fh.write("\t".join(s[col] for col in header) + "\n")
    print(f"Wrote {len(samples)} samples to {output}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Auto-generate samples.tsv from a FASTQ directory")
    parser.add_argument("--fastq-dir", required=True,
                        help="Directory containing FASTQ files")
    parser.add_argument("-o", "--output", default="config/samples.tsv",
                        help="Output sample sheet path")
    parser.add_argument("--groups", default="Sham,3day,7day,14day,28day",
                        help="Comma-separated group names")
    parser.add_argument("--reps", type=int, default=3,
                        help="Number of biological replicates per group")
    args = parser.parse_args()

    groups = [g.strip() for g in args.groups.split(",")]
    samples = find_fastqs(args.fastq_dir, groups, args.reps)

    if not samples:
        print("ERROR: No FASTQ files found. Check --fastq-dir and naming.",
              file=sys.stderr)
        sys.exit(1)

    write_tsv(samples, args.output)


if __name__ == "__main__":
    main()
