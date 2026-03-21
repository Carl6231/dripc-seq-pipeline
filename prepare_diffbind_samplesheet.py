#!/usr/bin/env python3
"""
prepare_diffbind_samplesheet.py
================================
Generate a DiffBind-compatible sample sheet CSV from the pipeline's samples.tsv.

DiffBind expects columns:
  SampleID, Tissue, Factor, Condition, Treatment, Replicate,
  bamReads, ControlID, bamControl, Peaks, PeakCaller

Usage:
  python scripts/prepare_diffbind_samplesheet.py \
      --samples config/samples.tsv \
      --peak-dir results/peaks \
      --bam-dir results/alignment \
      --output results/diffbind/diffbind_samplesheet.csv \
      --peak-caller stochhmm
"""

import argparse
import os
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Generate DiffBind sample sheet")
    parser.add_argument("--samples", required=True, help="Path to samples.tsv")
    parser.add_argument("--peak-dir", required=True, help="Directory with filtered BED peaks")
    parser.add_argument("--bam-dir", required=True, help="Directory with sorted BAMs")
    parser.add_argument("--output", required=True, help="Output CSV for DiffBind")
    parser.add_argument("--peak-caller", default="stochhmm",
                        choices=["stochhmm", "macs3"],
                        help="Peak caller used")
    args = parser.parse_args()

    samples = pd.read_csv(args.samples, sep="\t", comment="#", dtype=str)

    rows = []
    for _, row in samples.iterrows():
        sid = row["sample_id"]
        grp = row["group"]
        rep = row["rep"]
        has_input = row["has_input"].lower() in ("true", "yes", "1")

        bam_reads = os.path.abspath(
            os.path.join(args.bam_dir, f"{sid}.sorted.bam"))
        peaks = os.path.abspath(
            os.path.join(args.peak_dir, f"{sid}.filtered.bed"))

        bam_control = ""
        control_id = ""
        if has_input:
            bam_control = os.path.abspath(
                os.path.join(args.bam_dir, f"{sid}_input.sorted.bam"))
            control_id = f"{sid}_input"

        # DiffBind PeakCaller label
        peak_caller_label = "bed"

        rows.append({
            "SampleID": sid,
            "Tissue": "spinal_cord",
            "Factor": "R-loop",
            "Condition": grp,
            "Treatment": grp,
            "Replicate": rep,
            "bamReads": bam_reads,
            "ControlID": control_id,
            "bamControl": bam_control,
            "Peaks": peaks,
            "PeakCaller": peak_caller_label,
        })

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"DiffBind sample sheet written to {args.output} ({len(df)} samples)")


if __name__ == "__main__":
    main()
