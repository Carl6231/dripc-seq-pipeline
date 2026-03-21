#!/usr/bin/env bash
# =============================================================================
# run_bam_qc.sh — Generate samtools flagstat and stats for a BAM file
# =============================================================================
# Usage:
#   bash scripts/run_bam_qc.sh <input.bam> <output.flagstat> <output.stats>
# =============================================================================
set -euo pipefail

BAM="$1"
FLAGSTAT="$2"
STATS="$3"

if [ ! -f "$BAM" ]; then
    echo "ERROR: BAM file not found: $BAM" >&2
    exit 1
fi

echo "[$(date)] Generating flagstat for: $BAM"
samtools flagstat "$BAM" > "$FLAGSTAT"

echo "[$(date)] Generating stats for: $BAM"
samtools stats "$BAM" > "$STATS"

echo "[$(date)] BAM QC complete"
