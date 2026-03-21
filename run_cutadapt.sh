#!/usr/bin/env bash
# =============================================================================
# run_cutadapt.sh — Paired-end adapter/quality trimming with Cutadapt
# =============================================================================
# Usage:
#   bash scripts/run_cutadapt.sh \
#       <in_r1> <in_r2> <out_r1> <out_r2> \
#       <adapter_fwd> <adapter_rev> \
#       <quality_cutoff> <min_length> <threads> [extra_args]
# =============================================================================
set -euo pipefail

IN_R1="$1"
IN_R2="$2"
OUT_R1="$3"
OUT_R2="$4"
ADAPTER_FWD="$5"
ADAPTER_REV="$6"
QUALITY="$7"
MIN_LEN="$8"
THREADS="$9"
EXTRA="${10:-}"

if [ ! -f "$IN_R1" ]; then
    echo "ERROR: Input R1 not found: $IN_R1" >&2
    exit 1
fi
if [ ! -f "$IN_R2" ]; then
    echo "ERROR: Input R2 not found: $IN_R2" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUT_R1")"

echo "[$(date)] Cutadapt: trimming $IN_R1 + $IN_R2"

cutadapt \
    -a "$ADAPTER_FWD" \
    -A "$ADAPTER_REV" \
    -q "$QUALITY" \
    -m "$MIN_LEN" \
    -j "$THREADS" \
    $EXTRA \
    -o "$OUT_R1" \
    -p "$OUT_R2" \
    "$IN_R1" "$IN_R2"

echo "[$(date)] Cutadapt complete. Output: $OUT_R1, $OUT_R2"
