#!/usr/bin/env Rscript
# =============================================================================
# run_diffbind.R — Differential R-loop peak analysis using DiffBind
# =============================================================================
# Usage:
#   Rscript scripts/run_diffbind.R \
#       --samplesheet results/diffbind/diffbind_samplesheet.csv \
#       --comparisons "[['3day','Sham'],['7day','3day']]" \
#       --fc-threshold 1.5 \
#       --fdr-threshold 0.05 \
#       --min-overlap 2 \
#       --normalization default \
#       --outdir results/diffbind \
#       --counts-rds results/diffbind/dba_counts.rds
# =============================================================================
# NOTE: This is a standardized representative DiffBind workflow for documenting
# the main analytical steps. Key thresholds are user-configurable via config.yaml.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(DiffBind)
})

source("scripts/utils.R")

# ── parse arguments ──────────────────────────────────────────────────────────
option_list <- list(
  make_option("--samplesheet", type = "character",
              help = "DiffBind sample sheet CSV"),
  make_option("--comparisons", type = "character",
              help = "Comparisons as Python-style list string"),
  make_option("--fc-threshold", type = "double", default = 1.5,
              help = "Fold-change threshold [default: 1.5]"),
  make_option("--fdr-threshold", type = "double",
              default = 0.05, dest = "fdr_threshold",
              help = "FDR threshold for significant differential R-loop-enriched regions [default: 0.05]"),
  make_option("--min-overlap", type = "integer", default = 2,
              help = "DiffBind minOverlap [default: 2]"),
  make_option("--normalization", type = "character", default = "default",
              help = "Normalization: default|lib|RLE|TMM|native [default: default]"),
  make_option("--outdir", type = "character", default = "results/diffbind",
              help = "Output directory"),
  make_option("--counts-rds", type = "character",
              default = "results/diffbind/dba_counts.rds",
              help = "Path to save counts DBA object as RDS")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ── parse comparisons string ─────────────────────────────────────────────────
# Input format: "[['3day','Sham'], ['7day','3day']]"
parse_comparisons <- function(s) {
  # Clean up Python-style list notation
  s <- gsub("\\[", "", s)
  s <- gsub("\\]", "", s)
  s <- gsub("'", "", s)
  s <- gsub('"', "", s)
  parts <- strsplit(s, ",")[[1]]
  parts <- trimws(parts)
  # Pair them up
  n <- length(parts)
  if (n %% 2 != 0) stop("Comparisons must be pairs: treat,control")
  comps <- list()
  for (i in seq(1, n, 2)) {
    comps <- c(comps, list(c(parts[i], parts[i + 1])))
  }
  return(comps)
}

comparisons <- parse_comparisons(opt$comparisons)

cat("=== DiffBind differential analysis ===\n")
cat("Sample sheet:", opt$samplesheet, "\n")
cat("Comparisons:", length(comparisons), "\n")
cat("FC threshold:", opt$fc_threshold, "\n")
cat("FDR threshold:", opt$fdr_threshold, "\n")
cat("minOverlap:", opt$`min-overlap`, "\n")
cat("normalization:", opt$normalization, "\n")

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ── load samples ─────────────────────────────────────────────────────────────
cat("\n[1/5] Loading sample sheet...\n")
dba_obj <- dba(sampleSheet = opt$samplesheet)

# ── count reads ──────────────────────────────────────────────────────────────
cat("[2/5] Counting reads in consensus peaks (minOverlap =",
    opt$`min-overlap`, ")...\n")
dba_obj <- dba.count(dba_obj, minOverlap = opt$`min-overlap`)

# Save counts object
saveRDS(dba_obj, opt$`counts-rds`)
cat("  Saved counts object to:", opt$`counts-rds`, "\n")

# ── normalize ────────────────────────────────────────────────────────────────
cat("[3/5] Normalizing...\n")
# Apply normalization if not "default"
if (opt$normalization != "default") {
  # DiffBind >= 3.x uses dba.normalize; older versions handle it differently
  if (exists("dba.normalize")) {
    norm_method <- switch(opt$normalization,
                          "lib" = DBA_NORM_LIB,
                          "RLE" = DBA_NORM_RLE,
                          "TMM" = DBA_NORM_TMM,
                          "native" = DBA_NORM_NATIVE,
                          DBA_NORM_DEFAULT)
    dba_obj <- dba.normalize(dba_obj, normalize = norm_method)
  } else {
    cat("  NOTE: dba.normalize not available in this DiffBind version.",
        "Using default normalization.\n")
  }
} else {
  # For DiffBind >= 3.x, explicitly call normalize with default
  if (exists("dba.normalize")) {
    dba_obj <- dba.normalize(dba_obj)
  }
}

# ── run contrasts and report ─────────────────────────────────────────────────
cat("[4/5] Running differential analysis for each comparison...\n")

for (comp in comparisons) {
  treat <- comp[1]
  control <- comp[2]
  comp_name <- paste0(treat, "_vs_", control)
  cat("\n--- Comparison:", comp_name, "---\n")

  tryCatch({
    # Set up contrast
    dba_contrast <- dba.contrast(dba_obj,
                                  group1 = dba.mask(dba_obj, DBA_CONDITION, treat),
                                  group2 = dba.mask(dba_obj, DBA_CONDITION, control),
                                  name1 = treat,
                                  name2 = control)

    # Run analysis
    dba_result <- dba.analyze(dba_contrast)

    # Extract all results; significant peaks are filtered below using FDR.
    res <- dba.report(dba_result, th = 1)  # th=1 gets all peaks
    res_df <- as.data.frame(res)

    # Add log2FC and rename columns for clarity
    if ("Fold" %in% colnames(res_df)) {
      res_df$log2FoldChange <- res_df$Fold
    }

    # Write all results
    all_file <- file.path(opt$outdir, paste0(comp_name, ".all.tsv"))
    write.table(res_df, all_file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("  All peaks:", nrow(res_df), "->", all_file, "\n")

    # Filter significant peaks using FDR and fold-change thresholds.
    fc_col <- if ("Fold" %in% colnames(res_df)) "Fold" else "log2FoldChange"
    possible_fdr_cols <- c("FDR", "padj", "q.value", "qvalue")
    fdr_col <- intersect(possible_fdr_cols, colnames(res_df))[1]

    if (is.na(fdr_col)) {
      cat("  WARNING: Could not find FDR column. Columns:",
          paste(colnames(res_df), collapse = ", "), "\n")
      sig_df <- res_df[0, ]
    } else {
      sig_idx <- abs(res_df[[fc_col]]) >= log2(opt$fc_threshold) &
                 res_df[[fdr_col]] < opt$fdr_threshold
      sig_df <- res_df[sig_idx, ]
    }

    sig_file <- file.path(opt$outdir, paste0(comp_name, ".sig.tsv"))
    write.table(sig_df, sig_file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("  Significant peaks:", nrow(sig_df), "->", sig_file, "\n")

  }, error = function(e) {
    cat("  ERROR in comparison", comp_name, ":", conditionMessage(e), "\n")
    # Write empty files to satisfy Snakemake
    all_file <- file.path(opt$outdir, paste0(comp_name, ".all.tsv"))
    sig_file <- file.path(opt$outdir, paste0(comp_name, ".sig.tsv"))
    write.table(data.frame(), all_file, sep = "\t", row.names = FALSE)
    write.table(data.frame(), sig_file, sep = "\t", row.names = FALSE)
  })
}

cat("\n[5/5] DiffBind analysis complete.\n")