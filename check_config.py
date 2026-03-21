#!/usr/bin/env python3
"""
check_config.py
===============
Validate the config dictionary at Snakemake DAG-build time.
Raises informative errors for missing or suspicious settings.

Imported by the Snakefile:
    from check_config import validate_config
    validate_config(config)
"""

import sys


def validate_config(cfg):
    """Check required keys and warn about vendor-unknown parameters."""

    errors = []
    warnings = []

    # ── required top-level keys ───────────────────────────────────────────
    required_keys = ["samples_tsv", "genome", "cutadapt", "bowtie2",
                     "stochhmm", "macs3", "annotation", "diffbind",
                     "comparisons"]
    for k in required_keys:
        if k not in cfg:
            errors.append(f"Missing required config key: '{k}'")

    # ── genome ────────────────────────────────────────────────────────────
    genome = cfg.get("genome", {})
    if not genome.get("fasta"):
        warnings.append("genome.fasta is not set; some tools may need it")

    # ── Bowtie2 ───────────────────────────────────────────────────────────
    bt2 = cfg.get("bowtie2", {})
    if not bt2.get("index"):
        errors.append("bowtie2.index must be set to the Bowtie2 index prefix")

    # ── peak caller ───────────────────────────────────────────────────────
    peak_caller = cfg.get("peak_caller", "stochhmm")
    if peak_caller not in ("stochhmm", "macs3"):
        errors.append(f"peak_caller must be 'stochhmm' or 'macs3', got: {peak_caller}")

    # ── StochHMM ──────────────────────────────────────────────────────────
    if peak_caller == "stochhmm":
        shmm = cfg.get("stochhmm", {})
        if not shmm.get("binary"):
            errors.append(
                "stochhmm.binary is empty. The vendor used StochHMM v0.38 but did "
                "not disclose the exact binary path. Please set this to the path of "
                "your StochHMM executable, or switch peak_caller to 'macs3'."
            )
        if not shmm.get("model") and not shmm.get("command_template"):
            errors.append(
                "stochhmm.model is empty and no command_template is provided. "
                "The vendor did not disclose the HMM model file. Please provide "
                "either the model path or a full command_template."
            )

    # ── DiffBind — transparency about vendor-unknown params ───────────────
    db = cfg.get("diffbind", {})
    if db.get("normalization", "default") not in (
        "default", "lib", "RLE", "TMM", "native"
    ):
        errors.append(
            f"diffbind.normalization '{db['normalization']}' is not recognised. "
            "Accepted: default, lib, RLE, TMM, native"
        )
    # Warn about vendor-unknown params
    warnings.append(
        "NOTE: DiffBind minOverlap={}, normalization='{}' are user-configurable. "
        "The vendor (DiffBind 2.6.6) did not disclose these exact values.".format(
            db.get("minOverlap", 2), db.get("normalization", "default"))
    )

    # ── comparisons ───────────────────────────────────────────────────────
    comps = cfg.get("comparisons", {})
    if not comps.get("main"):
        errors.append("comparisons.main must contain at least one comparison")

    # ── annotation ────────────────────────────────────────────────────────
    ann = cfg.get("annotation", {})
    if not ann.get("txdb"):
        errors.append("annotation.txdb must be set")
    if not ann.get("orgdb"):
        errors.append("annotation.orgdb must be set")

    # ── report ────────────────────────────────────────────────────────────
    if errors:
        msg = "\n".join(f"  ERROR: {e}" for e in errors)
        sys.stderr.write(f"\n{'='*60}\nConfig validation FAILED:\n{msg}\n{'='*60}\n\n")
        raise ValueError("Config validation failed. See errors above.")

    if warnings:
        msg = "\n".join(f"  WARNING: {w}" for w in warnings)
        sys.stderr.write(f"\n{'='*60}\nConfig warnings:\n{msg}\n{'='*60}\n\n")
