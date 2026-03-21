#!/usr/bin/env bash
# =============================================================================
# run_bowtie2.sh — Paired-end alignment with Bowtie2 + sort/index with samtools
# =============================================================================
# Usage:
#   bash scripts/run_bowtie2.sh \
#       <r1.fq.gz> <r2.fq.gz> <output.sorted.bam> \
#       <bt2_index> <threads> [extra_args]
# =============================================================================
set -euo pipefail

R1="$1"
R2="$2"
OUT_BAM="$3"
INDEX="$4"
THREADS="$5"
EXTRA="${6:-}"

if [ ! -f "$R1" ]; then
    echo "ERROR: R1 not found: $R1" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUT_BAM")"

TMPBAM="${OUT_BAM%.sorted.bam}.tmp.bam"

echo "[$(date)] Bowtie2: aligning $R1 + $R2 to index $INDEX"

bowtie2 \
    -x "$INDEX" \
    -1 "$R1" \
    -2 "$R2" \
    --threads "$THREADS" \
    $EXTRA \
    2> >(tee "${OUT_BAM%.sorted.bam}.bowtie2.log" >&2) \
    | samtools view -@ "$THREADS" -bS -q 20 - \
    | samtools sort -@ "$THREADS" -o "$OUT_BAM" -

echo "[$(date)] Indexing BAM: $OUT_BAM"
samtools index "$OUT_BAM"

echo "[$(date)] Alignment complete: $OUT_BAM"
