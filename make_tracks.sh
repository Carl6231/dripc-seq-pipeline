#!/usr/bin/env bash
# =============================================================================
# make_tracks.sh — Generate RPM-normalized bigWig from BAM using deepTools
# =============================================================================
# Usage:
#   bash scripts/make_tracks.sh <input.bam> <output.bw> <norm> <binsize> <threads>
# =============================================================================
set -euo pipefail

BAM="$1"
BW="$2"
NORM="${3:-RPKM}"
BINSIZE="${4:-10}"
THREADS="${5:-4}"

if [ ! -f "$BAM" ]; then
    echo "ERROR: BAM file not found: $BAM" >&2
    exit 1
fi

mkdir -p "$(dirname "$BW")"

echo "[$(date)] bamCoverage: generating bigWig from $BAM"
echo "  Normalization: $NORM, binsize: $BINSIZE"

bamCoverage \
    -b "$BAM" \
    -o "$BW" \
    --normalizeUsing "$NORM" \
    --binSize "$BINSIZE" \
    -p "$THREADS" \
    --extendReads

echo "[$(date)] bigWig track created: $BW"
