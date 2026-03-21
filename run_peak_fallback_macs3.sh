#!/usr/bin/env bash
# =============================================================================
# run_peak_fallback_macs3.sh — MACS3 broad peak calling (fallback)
# =============================================================================
# Usage:
#   bash scripts/run_peak_fallback_macs3.sh \
#       <ip_bam> <input_bam_or_empty> <output_bed> \
#       <genome_size> <broad> <qvalue> <extra_args> <sample_id>
# =============================================================================
set -euo pipefail

IP_BAM="$1"
INPUT_BAM="$2"
OUTPUT_BED="$3"
GSIZE="$4"
BROAD="$5"
QVALUE="$6"
EXTRA="$7"
SAMPLE_ID="$8"

if [ ! -f "$IP_BAM" ]; then
    echo "ERROR: IP BAM not found: $IP_BAM" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT_BED")"

TMPDIR="results/peaks/macs3_tmp_${SAMPLE_ID}"
mkdir -p "$TMPDIR"

echo "[$(date)] MACS3 peak calling for $SAMPLE_ID"

# Build MACS3 command
CMD="macs3 callpeak"
CMD="$CMD -t $IP_BAM"
CMD="$CMD -f BAMPE"
CMD="$CMD -g $GSIZE"
CMD="$CMD -q $QVALUE"
CMD="$CMD -n $SAMPLE_ID"
CMD="$CMD --outdir $TMPDIR"

# Add input control if available
if [ -n "$INPUT_BAM" ] && [ -f "$INPUT_BAM" ]; then
    CMD="$CMD -c $INPUT_BAM"
    echo "  Using input control: $INPUT_BAM"
else
    echo "  No input control — running without control BAM"
fi

# Broad peak mode
if [ "$BROAD" = "true" ] || [ "$BROAD" = "True" ]; then
    CMD="$CMD --broad"
    PEAK_SUFFIX="broadPeak"
else
    PEAK_SUFFIX="narrowPeak"
fi

# Extra arguments
if [ -n "$EXTRA" ]; then
    CMD="$CMD $EXTRA"
fi

echo "  Command: $CMD"
eval "$CMD"

# Convert peak file to simple BED
PEAK_FILE="${TMPDIR}/${SAMPLE_ID}_peaks.${PEAK_SUFFIX}"
if [ ! -f "$PEAK_FILE" ]; then
    echo "WARNING: MACS3 peak file not found: $PEAK_FILE" >&2
    # Try the other suffix
    if [ "$PEAK_SUFFIX" = "broadPeak" ]; then
        PEAK_FILE="${TMPDIR}/${SAMPLE_ID}_peaks.narrowPeak"
    else
        PEAK_FILE="${TMPDIR}/${SAMPLE_ID}_peaks.broadPeak"
    fi
fi

if [ -f "$PEAK_FILE" ]; then
    # Extract chrom, start, end, name, score, strand
    cut -f1-6 "$PEAK_FILE" | sort -k1,1 -k2,2n > "$OUTPUT_BED"
else
    echo "WARNING: No MACS3 peak file found. Creating empty output." >&2
    touch "$OUTPUT_BED"
fi

# Clean up temp dir
rm -rf "$TMPDIR"

NPEAKS=$(wc -l < "$OUTPUT_BED")
echo "[$(date)] MACS3 called $NPEAKS raw peaks for $SAMPLE_ID"
echo "[$(date)] Output: $OUTPUT_BED"
