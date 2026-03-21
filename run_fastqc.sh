#!/usr/bin/env bash
# =============================================================================
# run_fastqc.sh — Run FastQC on a single FASTQ file
# =============================================================================
# Usage:
#   bash scripts/run_fastqc.sh <fastq> <outdir> <threads>
# =============================================================================
set -euo pipefail

FASTQ="$1"
OUTDIR="$2"
THREADS="${3:-2}"

if [ ! -f "$FASTQ" ]; then
    echo "ERROR: FASTQ file not found: $FASTQ" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

echo "[$(date)] Running FastQC on: $FASTQ"
fastqc -t "$THREADS" -o "$OUTDIR" "$FASTQ"
echo "[$(date)] FastQC complete for: $FASTQ"
