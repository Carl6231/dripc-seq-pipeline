"""
Snakefile — DRIPc-seq analysis pipeline for mouse spinal cord injury (SCI)
==========================================================================
Groups : Sham, 3day, 7day, 14day, 28day  (3 biological replicates each)
Platform: Illumina NovaSeq 6000, PE 2x150
Genome  : GRCm38

Main flow:
  FastQC → Cutadapt → Bowtie2 → StochHMM (or MACS3 fallback)
  → filter peaks ≥100 bp → DiffBind → annotation → enrichment

Usage:
  snakemake --use-conda -j 16 --configfile config/config.yaml
"""

import pandas as pd
import os
import sys

# ── load & validate config ────────────────────────────────────────────────────
configfile: "config/config.yaml"

# run config checker at DAG build time
sys.path.insert(0, "scripts")
from check_config import validate_config
validate_config(config)

# ── read sample sheet ─────────────────────────────────────────────────────────
samples_tsv = config["samples_tsv"]
if not os.path.exists(samples_tsv):
    raise FileNotFoundError(
        f"Sample sheet not found: {samples_tsv}\n"
        "Run: python scripts/build_samplesheet.py --fastq-dir <dir> -o config/samples.tsv"
    )

SAMPLES = pd.read_csv(samples_tsv, sep="\t", comment="#", dtype=str).set_index(
    "sample_id", drop=False
)
SAMPLE_IDS = SAMPLES["sample_id"].tolist()
GROUPS = sorted(SAMPLES["group"].unique())

# check which samples have input controls
HAS_INPUT = {sid: row["has_input"].lower() in ("true", "yes", "1")
             for sid, row in SAMPLES.iterrows()}

# ── peak caller selection ─────────────────────────────────────────────────────
PEAK_CALLER = config.get("peak_caller", "stochhmm")  # "stochhmm" or "macs3"

# ── comparison design ─────────────────────────────────────────────────────────
MAIN_COMPARISONS = config["comparisons"]["main"]      # list of [treat, control]
OPTIONAL_COMPARISONS = config["comparisons"].get("optional", [])
if config.get("run_optional_comparisons", False):
    ALL_COMPARISONS = MAIN_COMPARISONS + OPTIONAL_COMPARISONS
else:
    ALL_COMPARISONS = MAIN_COMPARISONS

COMP_NAMES = [f"{c[0]}_vs_{c[1]}" for c in ALL_COMPARISONS]

# ── DEG overlap ───────────────────────────────────────────────────────────────
DEG_FILE = config.get("deg_file", "")
RUN_DEG_OVERLAP = bool(DEG_FILE) and os.path.isfile(DEG_FILE)

# ── helper: list all FASTQ files to QC ────────────────────────────────────────
def all_fastq_for_qc():
    """Return list of FASTQ paths (IP + input if present) for FastQC."""
    fqs = []
    for _, row in SAMPLES.iterrows():
        fqs += [row["ip_fq1"], row["ip_fq2"]]
        if row["has_input"].lower() in ("true", "yes", "1"):
            fqs += [row["input_fq1"], row["input_fq2"]]
    return [f for f in fqs if pd.notna(f) and f != "" and f != "NA"]

ALL_FASTQS = all_fastq_for_qc()

# ── targets ───────────────────────────────────────────────────────────────────
rule all:
    input:
        # QC
        "results/qc/qc_summary.tsv",
        # Alignment QC
        expand("results/alignment/{sid}.flagstat", sid=SAMPLE_IDS),
        # Tracks
        expand("results/tracks/{sid}.bw", sid=SAMPLE_IDS),
        # Filtered peaks
        expand("results/peaks/{sid}.filtered.bed", sid=SAMPLE_IDS),
        # Annotation
        expand("results/annotation/{sid}.annotated.tsv", sid=SAMPLE_IDS),
        # DiffBind
        expand("results/diffbind/{comp}.all.tsv", comp=COMP_NAMES),
        expand("results/plots/{comp}.volcano.pdf", comp=COMP_NAMES),
        # Heatmaps
        "results/plots/sample_correlation_heatmap.pdf",
        "results/plots/diff_peaks_heatmap.pdf",
        # Enrichment
        expand("results/enrichment/{comp}.go_bp.tsv", comp=COMP_NAMES),
        # DEG overlap (conditional)
        (expand("results/overlap/{comp}.overlap_summary.tsv", comp=COMP_NAMES)
         if RUN_DEG_OVERLAP else []),


# ══════════════════════════════════════════════════════════════════════════════
# RULE: FastQC
# ══════════════════════════════════════════════════════════════════════════════
rule fastqc:
    input:
        fq=lambda wc: wc.fqpath,
    output:
        html="results/qc/fastqc/{fqbase}_fastqc.html",
        zipf="results/qc/fastqc/{fqbase}_fastqc.zip",
    params:
        outdir="results/qc/fastqc",
        fqpath=lambda wc: wc.fqpath,
    log:
        "logs/fastqc/{fqbase}.log",
    conda:
        "envs/base.yaml"
    threads: 2
    shell:
        """
        bash scripts/run_fastqc.sh \
            {params.fqpath} {params.outdir} {threads} \
            > {log} 2>&1
        """

# Because Snakemake cannot easily wildcard over the full FASTQ paths,
# we instead use a checkpoint / aggregate approach via the summary rule.

rule fastqc_all:
    """Run FastQC on every FASTQ file and then produce a summary."""
    input:
        fqs=ALL_FASTQS,
    output:
        summary="results/qc/qc_summary.tsv",
    params:
        outdir="results/qc/fastqc",
    log:
        "logs/fastqc/fastqc_all.log",
    conda:
        "envs/base.yaml"
    threads: config.get("fastqc_threads", 4)
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.fqs} > {log} 2>&1

        # Build summary
        echo -e "file\\ttotal_sequences\\tpoor_quality\\tsequence_length\\tpct_gc" \
            > {output.summary}
        for z in {params.outdir}/*_fastqc.zip; do
            base=$(basename "$z" _fastqc.zip)
            unzip -p "$z" "$base"_fastqc/fastqc_data.txt 2>/dev/null | awk '
                BEGIN{{OFS="\\t"; ts=""; pq=""; sl=""; gc=""}}
                /^Total Sequences/{{ts=$2}}
                /^Sequences flagged as poor quality/{{pq=$NF}}
                /^Sequence length/{{sl=$NF}}
                /^%GC/{{gc=$2}}
                END{{print ENVIRON["base"], ts, pq, sl, gc}}
            ' base="$base" >> {output.summary} || true
        done
        echo "[$(date)] FastQC summary written to {output.summary}" >> {log}
        """

# ══════════════════════════════════════════════════════════════════════════════
# RULE: Cutadapt
# ══════════════════════════════════════════════════════════════════════════════
rule cutadapt:
    input:
        r1=lambda wc: SAMPLES.loc[wc.sid, "ip_fq1"],
        r2=lambda wc: SAMPLES.loc[wc.sid, "ip_fq2"],
    output:
        r1="results/trimmed/{sid}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sid}_R2.trimmed.fastq.gz",
    params:
        adapter_fwd=config["cutadapt"]["adapter_fwd"],
        adapter_rev=config["cutadapt"]["adapter_rev"],
        quality=config["cutadapt"]["quality_cutoff"],
        min_len=config["cutadapt"]["min_length"],
        extra=config["cutadapt"].get("extra_args", ""),
    log:
        "logs/cutadapt/{sid}.log",
    conda:
        "envs/base.yaml"
    threads: config.get("cutadapt_threads", 4)
    shell:
        """
        bash scripts/run_cutadapt.sh \
            {input.r1} {input.r2} \
            {output.r1} {output.r2} \
            "{params.adapter_fwd}" "{params.adapter_rev}" \
            {params.quality} {params.min_len} \
            {threads} "{params.extra}" \
            > {log} 2>&1
        """

rule cutadapt_input:
    """Trim input-control FASTQ (only when present)."""
    input:
        r1=lambda wc: SAMPLES.loc[wc.sid, "input_fq1"],
        r2=lambda wc: SAMPLES.loc[wc.sid, "input_fq2"],
    output:
        r1="results/trimmed/{sid}_input_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sid}_input_R2.trimmed.fastq.gz",
    params:
        adapter_fwd=config["cutadapt"]["adapter_fwd"],
        adapter_rev=config["cutadapt"]["adapter_rev"],
        quality=config["cutadapt"]["quality_cutoff"],
        min_len=config["cutadapt"]["min_length"],
        extra=config["cutadapt"].get("extra_args", ""),
    log:
        "logs/cutadapt/{sid}_input.log",
    conda:
        "envs/base.yaml"
    threads: config.get("cutadapt_threads", 4)
    shell:
        """
        bash scripts/run_cutadapt.sh \
            {input.r1} {input.r2} \
            {output.r1} {output.r2} \
            "{params.adapter_fwd}" "{params.adapter_rev}" \
            {params.quality} {params.min_len} \
            {threads} "{params.extra}" \
            > {log} 2>&1
        """

# ══════════════════════════════════════════════════════════════════════════════
# RULE: Bowtie2 alignment
# ══════════════════════════════════════════════════════════════════════════════
rule bowtie2_ip:
    input:
        r1="results/trimmed/{sid}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sid}_R2.trimmed.fastq.gz",
    output:
        bam="results/alignment/{sid}.sorted.bam",
        bai="results/alignment/{sid}.sorted.bam.bai",
    params:
        index=config["bowtie2"]["index"],
        extra=config["bowtie2"].get("extra_args", ""),
    log:
        "logs/bowtie2/{sid}.log",
    conda:
        "envs/base.yaml"
    threads: config.get("bowtie2_threads", 8)
    shell:
        """
        bash scripts/run_bowtie2.sh \
            {input.r1} {input.r2} \
            {output.bam} \
            {params.index} \
            {threads} "{params.extra}" \
            > {log} 2>&1
        """

rule bowtie2_input:
    """Align input-control reads (only when has_input is true)."""
    input:
        r1="results/trimmed/{sid}_input_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sid}_input_R2.trimmed.fastq.gz",
    output:
        bam="results/alignment/{sid}_input.sorted.bam",
        bai="results/alignment/{sid}_input.sorted.bam.bai",
    params:
        index=config["bowtie2"]["index"],
        extra=config["bowtie2"].get("extra_args", ""),
    log:
        "logs/bowtie2/{sid}_input.log",
    conda:
        "envs/base.yaml"
    threads: config.get("bowtie2_threads", 8)
    shell:
        """
        bash scripts/run_bowtie2.sh \
            {input.r1} {input.r2} \
            {output.bam} \
            {params.index} \
            {threads} "{params.extra}" \
            > {log} 2>&1
        """

# ══════════════════════════════════════════════════════════════════════════════
# RULE: BAM QC (flagstat + stats)
# ══════════════════════════════════════════════════════════════════════════════
rule bam_qc:
    input:
        bam="results/alignment/{sid}.sorted.bam",
    output:
        flagstat="results/alignment/{sid}.flagstat",
        stats="results/alignment/{sid}.stats",
    log:
        "logs/bam_qc/{sid}.log",
    conda:
        "envs/base.yaml"
    shell:
        """
        bash scripts/run_bam_qc.sh {input.bam} \
            {output.flagstat} {output.stats} \
            > {log} 2>&1
        """

# ══════════════════════════════════════════════════════════════════════════════
# RULE: signal tracks (bigWig)
# ══════════════════════════════════════════════════════════════════════════════
rule make_tracks:
    input:
        bam="results/alignment/{sid}.sorted.bam",
        bai="results/alignment/{sid}.sorted.bam.bai",
    output:
        bw="results/tracks/{sid}.bw",
    params:
        norm=config.get("track_normalization", "RPKM"),
        binsize=config.get("track_binsize", 10),
    log:
        "logs/tracks/{sid}.log",
    conda:
        "envs/base.yaml"
    threads: config.get("track_threads", 4)
    shell:
        """
        bash scripts/make_tracks.sh \
            {input.bam} {output.bw} \
            {params.norm} {params.binsize} {threads} \
            > {log} 2>&1
        """

# ══════════════════════════════════════════════════════════════════════════════
# RULE: Peak calling — StochHMM (primary)
# ══════════════════════════════════════════════════════════════════════════════
rule stochhmm_peaks:
    input:
        bam="results/alignment/{sid}.sorted.bam",
        # conditionally include input BAM
        input_bam=lambda wc: (
            "results/alignment/{}_input.sorted.bam".format(wc.sid)
            if HAS_INPUT.get(wc.sid, False) else []
        ),
    output:
        bed="results/peaks/{sid}.stochhmm.raw.bed",
    params:
        binary=config["stochhmm"]["binary"],
        model=config["stochhmm"]["model"],
        extra=config["stochhmm"].get("extra_args", ""),
        template=config["stochhmm"].get("command_template", ""),
        has_input=lambda wc: "true" if HAS_INPUT.get(wc.sid, False) else "false",
    log:
        "logs/peaks/{sid}.stochhmm.log",
    conda:
        "envs/base.yaml"
    shell:
        """
        input_bam=""
        if [ "{params.has_input}" = "true" ]; then
            input_bam="results/alignment/{wildcards.sid}_input.sorted.bam"
        fi

        bash scripts/run_stochhmm_wrapper.sh \
            {input.bam} \
            "$input_bam" \
            {output.bed} \
            "{params.binary}" \
            "{params.model}" \
            "{params.extra}" \
            "{params.template}" \
            > {log} 2>&1
        """

# ══════════════════════════════════════════════════════════════════════════════
# RULE: Peak calling — MACS3 fallback
# ══════════════════════════════════════════════════════════════════════════════
rule macs3_peaks:
    input:
        bam="results/alignment/{sid}.sorted.bam",
        input_bam=lambda wc: (
            "results/alignment/{}_input.sorted.bam".format(wc.sid)
            if HAS_INPUT.get(wc.sid, False) else []
        ),
    output:
        bed="results/peaks/{sid}.macs3.raw.bed",
    params:
        genome_size=config["macs3"].get("genome_size", "mm"),
        broad=config["macs3"].get("broad", True),
        qvalue=config["macs3"].get("qvalue", 0.05),
        extra=config["macs3"].get("extra_args", ""),
        has_input=lambda wc: "true" if HAS_INPUT.get(wc.sid, False) else "false",
    log:
        "logs/peaks/{sid}.macs3.log",
    conda:
        "envs/macs3.yaml"
    shell:
        """
        input_bam=""
        if [ "{params.has_input}" = "true" ]; then
            input_bam="results/alignment/{wildcards.sid}_input.sorted.bam"
        fi

        bash scripts/run_peak_fallback_macs3.sh \
            {input.bam} \
            "$input_bam" \
            {output.bed} \
            "{params.genome_size}" \
            "{params.broad}" \
            "{params.qvalue}" \
            "{params.extra}" \
            "{wildcards.sid}" \
            > {log} 2>&1
        """

# ══════════════════════════════════════════════════════════════════════════════
# RULE: filter peaks by width (≥100 bp)
# ══════════════════════════════════════════════════════════════════════════════
def raw_peak_file(wc):
    """Return the raw peak file depending on peak_caller selection."""
    if PEAK_CALLER == "macs3":
        return f"results/peaks/{wc.sid}.macs3.raw.bed"
    return f"results/peaks/{wc.sid}.stochhmm.raw.bed"

rule filter_peaks:
    input:
        bed=raw_peak_file,
    output:
        bed="results/peaks/{sid}.filtered.bed",
    params:
        min_width=config.get("peak_min_width", 100),
    log:
        "logs/peaks/{sid}.filter.log",
    conda:
        "envs/base.yaml"
    shell:
        """
        python scripts/filter_peaks_by_width.py \
            --input {input.bed} \
            --output {output.bed} \
            --min-width {params.min_width} \
            > {log} 2>&1
        """

# ══════════════════════════════════════════════════════════════════════════════
# RULE: annotate peaks
# ══════════════════════════════════════════════════════════════════════════════
rule annotate_peaks:
    input:
        bed="results/peaks/{sid}.filtered.bed",
    output:
        tsv="results/annotation/{sid}.annotated.tsv",
    params:
        txdb=config["annotation"]["txdb"],
        orgdb=config["annotation"]["orgdb"],
        promoter_up=config["annotation"].get("promoter_upstream", 2000),
        promoter_down=config["annotation"].get("promoter_downstream", 200),
        terminator_down=config["annotation"].get("terminator_downstream", 2000),
    log:
        "logs/annotation/{sid}.log",
    conda:
        "envs/r.yaml"
    shell:
        """
        Rscript scripts/annotate_peaks.R \
            --peaks {input.bed} \
            --output {output.tsv} \
            --txdb {params.txdb} \
            --orgdb {params.orgdb} \
            --promoter-up {params.promoter_up} \
            --promoter-down {params.promoter_down} \
            --terminator-down {params.terminator_down} \
            > {log} 2>&1
        """

# ══════════════════════════════════════════════════════════════════════════════
# RULE: prepare DiffBind sample sheet
# ══════════════════════════════════════════════════════════════════════════════
rule prepare_diffbind_sheet:
    input:
        peaks=expand("results/peaks/{sid}.filtered.bed", sid=SAMPLE_IDS),
        bams=expand("results/alignment/{sid}.sorted.bam", sid=SAMPLE_IDS),
    output:
        sheet="results/diffbind/diffbind_samplesheet.csv",
    params:
        samples_tsv=samples_tsv,
        peak_caller=PEAK_CALLER,
    log:
        "logs/diffbind/prepare_sheet.log",
    conda:
        "envs/base.yaml"
    shell:
        """
        python scripts/prepare_diffbind_samplesheet.py \
            --samples {params.samples_tsv} \
            --peak-dir results/peaks \
            --bam-dir results/alignment \
            --output {output.sheet} \
            --peak-caller {params.peak_caller} \
            > {log} 2>&1
        """

# ══════════════════════════════════════════════════════════════════════════════
# RULE: DiffBind
# ══════════════════════════════════════════════════════════════════════════════
rule diffbind:
    input:
        sheet="results/diffbind/diffbind_samplesheet.csv",
    output:
        all_tsv=expand("results/diffbind/{comp}.all.tsv", comp=COMP_NAMES),
        sig_tsv=expand("results/diffbind/{comp}.sig.tsv", comp=COMP_NAMES),
        counts="results/diffbind/dba_counts.rds",
    params:
        comparisons=str(ALL_COMPARISONS),
        fc_threshold=config["diffbind"].get("fc_threshold", 1.5),
        fdr_threshold=config["diffbind"].get("fdr_threshold", 0.05),
        min_overlap=config["diffbind"].get("minOverlap", 2),
        normalization=config["diffbind"].get("normalization", "default"),
        outdir="results/diffbind",
    log:
        "logs/diffbind/diffbind.log",
    conda:
        "envs/r.yaml"
    shell:
        """
        Rscript scripts/run_diffbind.R \
            --samplesheet {input.sheet} \
            --comparisons '{params.comparisons}' \
            --fc-threshold {params.fc_threshold} \
            --fdr-threshold {params.fdr_threshold} \
            --min-overlap {params.min_overlap} \
            --normalization "{params.normalization}" \
            --outdir {params.outdir} \
            --counts-rds {output.counts} \
            > {log} 2>&1
        """

# ══════════════════════════════════════════════════════════════════════════════
# RULE: volcano plots
# ══════════════════════════════════════════════════════════════════════════════
rule volcano:
    input:
        tsv="results/diffbind/{comp}.all.tsv",
    output:
        pdf="results/plots/{comp}.volcano.pdf",
        png="results/plots/{comp}.volcano.png",
    params:
        fc_threshold=config["diffbind"].get("fc_threshold", 1.5),
        fdr_threshold=config["diffbind"].get("fdr_threshold", 0.05),
    log:
        "logs/plots/{comp}.volcano.log",
    conda:
        "envs/r.yaml"
    shell:
        """
        Rscript scripts/plot_volcano.R \
            --input {input.tsv} \
            --output-pdf {output.pdf} \
            --output-png {output.png} \
            --fc-threshold {params.fc_threshold} \
            --fdr-threshold {params.fdr_threshold} \
            --title "{wildcards.comp}" \
            > {log} 2>&1
        """

# ══════════════════════════════════════════════════════════════════════════════
# RULE: heatmaps
# ══════════════════════════════════════════════════════════════════════════════
rule heatmap:
    input:
        counts="results/diffbind/dba_counts.rds",
        sig_files=expand("results/diffbind/{comp}.sig.tsv", comp=COMP_NAMES),
    output:
        corr_pdf="results/plots/sample_correlation_heatmap.pdf",
        corr_png="results/plots/sample_correlation_heatmap.png",
        diff_pdf="results/plots/diff_peaks_heatmap.pdf",
        diff_png="results/plots/diff_peaks_heatmap.png",
    params:
        sig_dir="results/diffbind",
    log:
        "logs/plots/heatmap.log",
    conda:
        "envs/r.yaml"
    shell:
        """
        Rscript scripts/plot_heatmap.R \
            --counts-rds {input.counts} \
            --sig-dir {params.sig_dir} \
            --output-corr-pdf {output.corr_pdf} \
            --output-corr-png {output.corr_png} \
            --output-diff-pdf {output.diff_pdf} \
            --output-diff-png {output.diff_png} \
            > {log} 2>&1
        """

# ══════════════════════════════════════════════════════════════════════════════
# RULE: DEG overlap (conditional)
# ══════════════════════════════════════════════════════════════════════════════
if RUN_DEG_OVERLAP:
    rule deg_overlap:
        input:
            sig="results/diffbind/{comp}.sig.tsv",
            deg=DEG_FILE,
        output:
            summary="results/overlap/{comp}.overlap_summary.tsv",
            genes="results/overlap/{comp}.overlap_genes.tsv",
            venn="results/overlap/{comp}.venn.pdf",
        log:
            "logs/overlap/{comp}.log",
        conda:
            "envs/r.yaml"
        shell:
            """
            Rscript scripts/deg_overlap.R \
                --sig-peaks {input.sig} \
                --deg-file {input.deg} \
                --output-summary {output.summary} \
                --output-genes {output.genes} \
                --output-venn {output.venn} \
                --comparison "{wildcards.comp}" \
                > {log} 2>&1
            """

# ══════════════════════════════════════════════════════════════════════════════
# RULE: enrichment (GO / KEGG)
# ══════════════════════════════════════════════════════════════════════════════
rule enrichment:
    input:
        sig="results/diffbind/{comp}.sig.tsv",
    output:
        go_bp="results/enrichment/{comp}.go_bp.tsv",
        go_mf="results/enrichment/{comp}.go_mf.tsv",
        kegg="results/enrichment/{comp}.kegg.tsv",
        go_pdf="results/enrichment/{comp}.go_bp.pdf",
        kegg_pdf="results/enrichment/{comp}.kegg.pdf",
    params:
        orgdb=config["annotation"]["orgdb"],
        pval_cutoff=config.get("enrichment_pval", 0.05),
        min_genes=config.get("enrichment_min_genes", 5),
    log:
        "logs/enrichment/{comp}.log",
    conda:
        "envs/r.yaml"
    shell:
        """
        Rscript scripts/run_enrichment.R \
            --sig-peaks {input.sig} \
            --orgdb {params.orgdb} \
            --pval-cutoff {params.pval_cutoff} \
            --min-genes {params.min_genes} \
            --output-go-bp {output.go_bp} \
            --output-go-mf {output.go_mf} \
            --output-kegg {output.kegg} \
            --output-go-pdf {output.go_pdf} \
            --output-kegg-pdf {output.kegg_pdf} \
            --comparison "{wildcards.comp}" \
            > {log} 2>&1
        """
