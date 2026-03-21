#!/usr/bin/env bash
# =============================================================================
# run_stochhmm_wrapper.sh — Wrapper for StochHMM peak calling
# =============================================================================
# IMPORTANT: The vendor used StochHMM v0.38 for DRIPc-seq peak calling but
# did NOT disclose the exact HMM model, emission parameters, or CLI syntax.
#
# This wrapper requires the user to configure:
#   stochhmm.binary           — path to the StochHMM executable
#   stochhmm.model            — path to the HMM model file
#   stochhmm.extra_args       — any additional arguments
#   stochhmm.command_template — (optional) full command template
#
# If a command_template is provided, it is used as-is with variable
# substitution. Otherwise, a default invocation is constructed.
# =============================================================================
# Usage:
#   bash scripts/run_stochhmm_wrapper.sh \
#       <ip_bam> <input_bam_or_empty> <output_bed> \
#       <binary> <model> <extra_args> <command_template>
# =============================================================================
set -euo pipefail

IP_BAM="$1"
INPUT_BAM="$2"          # may be empty string
OUTPUT_BED="$3"
BINARY="$4"
MODEL="$5"
EXTRA_ARGS="$6"
CMD_TEMPLATE="$7"

# ── validate required parameters ─────────────────────────────────────────────
if [ -z "$BINARY" ]; then
    echo "============================================================" >&2
    echo "ERROR: stochhmm.binary is not configured." >&2
    echo "" >&2
    echo "The vendor used StochHMM v0.38 but did not disclose the" >&2
    echo "exact binary path or installation. Please either:" >&2
    echo "  1) Set stochhmm.binary in config/config.yaml" >&2
    echo "  2) Switch to peak_caller: 'macs3' in config/config.yaml" >&2
    echo "============================================================" >&2
    exit 1
fi

if [ ! -x "$BINARY" ] && ! command -v "$BINARY" &>/dev/null; then
    echo "ERROR: StochHMM binary not found or not executable: $BINARY" >&2
    exit 1
fi

if [ -z "$MODEL" ] && [ -z "$CMD_TEMPLATE" ]; then
    echo "============================================================" >&2
    echo "ERROR: stochhmm.model is not configured and no" >&2
    echo "command_template is provided." >&2
    echo "" >&2
    echo "The vendor did not disclose the HMM model parameters." >&2
    echo "Please provide either:" >&2
    echo "  1) stochhmm.model — path to the .hmm model file" >&2
    echo "  2) stochhmm.command_template — full command with placeholders" >&2
    echo "============================================================" >&2
    exit 1
fi

if [ -n "$MODEL" ] && [ ! -f "$MODEL" ]; then
    echo "ERROR: StochHMM model file not found: $MODEL" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT_BED")"

echo "[$(date)] StochHMM peak calling"
echo "  IP BAM:    $IP_BAM"
echo "  Input BAM: ${INPUT_BAM:-none}"
echo "  Binary:    $BINARY"
echo "  Model:     ${MODEL:-'(using command_template)'}"

# ── build and execute command ────────────────────────────────────────────────
if [ -n "$CMD_TEMPLATE" ]; then
    # Use user-provided command template
    # Available placeholders: {binary} {model} {bam} {input_bam} {output} {extra}
    CMD="$CMD_TEMPLATE"
    CMD="${CMD//\{binary\}/$BINARY}"
    CMD="${CMD//\{model\}/$MODEL}"
    CMD="${CMD//\{bam\}/$IP_BAM}"
    CMD="${CMD//\{input_bam\}/$INPUT_BAM}"
    CMD="${CMD//\{output\}/$OUTPUT_BED}"
    CMD="${CMD//\{extra\}/$EXTRA_ARGS}"
    echo "[$(date)] Executing command template:"
    echo "  $CMD"
    eval "$CMD"
else
    # Construct default command
    # NOTE: This is a best-guess invocation. The exact StochHMM CLI for
    # DRIPc-seq peak calling was not disclosed by the vendor. Users should
    # verify and adjust the command or provide a command_template.
    CMD="$BINARY"
    CMD="$CMD -model $MODEL"
    CMD="$CMD -seq $IP_BAM"
    if [ -n "$INPUT_BAM" ] && [ -f "$INPUT_BAM" ]; then
        CMD="$CMD -input $INPUT_BAM"
    fi
    CMD="$CMD $EXTRA_ARGS"

    echo "[$(date)] Executing constructed command:"
    echo "  $CMD"
    echo "  NOTE: This is a best-guess invocation. The vendor did not"
    echo "  disclose exact StochHMM CLI parameters. Please verify output."

    # Run StochHMM and capture raw output
    TMPOUT="${OUTPUT_BED}.raw.tmp"
    eval "$CMD" > "$TMPOUT"

    # Convert StochHMM output to BED format
    # NOTE: The exact output format of the vendor's StochHMM run is unknown.
    # Below assumes a tab-delimited output with chrom, start, end (at minimum).
    # If your StochHMM outputs a different format, you may need to adjust this.
    if [ -f "$TMPOUT" ]; then
        # Attempt to extract BED-like columns (chrom, start, end)
        awk 'BEGIN{OFS="\t"} NF>=3 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ {
            print $1, $2, $3
        }' "$TMPOUT" | sort -k1,1 -k2,2n > "$OUTPUT_BED"
        rm -f "$TMPOUT"
    else
        echo "ERROR: StochHMM produced no output" >&2
        exit 1
    fi
fi

# ── verify output ────────────────────────────────────────────────────────────
if [ ! -s "$OUTPUT_BED" ]; then
    echo "WARNING: StochHMM output is empty: $OUTPUT_BED" >&2
    # Create empty BED to avoid downstream errors
    touch "$OUTPUT_BED"
fi

NPEAKS=$(wc -l < "$OUTPUT_BED")
echo "[$(date)] StochHMM called $NPEAKS raw peaks"
echo "[$(date)] Output: $OUTPUT_BED"
