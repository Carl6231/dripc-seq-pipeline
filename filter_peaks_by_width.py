#!/usr/bin/env python3
"""
filter_peaks_by_width.py
========================
Remove peaks shorter than a specified minimum width (default: 100 bp).

Usage:
  python scripts/filter_peaks_by_width.py \
      --input peaks.raw.bed \
      --output peaks.filtered.bed \
      --min-width 100
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Filter BED peaks by minimum width")
    parser.add_argument("--input", required=True, help="Input BED file")
    parser.add_argument("--output", required=True, help="Output BED file")
    parser.add_argument("--min-width", type=int, default=100,
                        help="Minimum peak width in bp (default: 100)")
    args = parser.parse_args()

    total = 0
    kept = 0
    removed = 0

    with open(args.input) as fin, open(args.output, "w") as fout:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                fout.write(line + "\n")
                continue
            total += 1
            fields = line.split("\t")
            if len(fields) < 3:
                print(f"WARNING: Skipping malformed line: {line}",
                      file=sys.stderr)
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                print(f"WARNING: Non-integer coordinates: {line}",
                      file=sys.stderr)
                continue

            width = end - start
            if width >= args.min_width:
                fout.write(line + "\n")
                kept += 1
            else:
                removed += 1

    print(f"[filter_peaks] Total: {total}, Kept (>={args.min_width} bp): {kept}, "
          f"Removed (<{args.min_width} bp): {removed}")


if __name__ == "__main__":
    main()
