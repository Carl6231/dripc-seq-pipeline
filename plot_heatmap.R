#!/usr/bin/env Rscript
# =============================================================================
# plot_heatmap.R — Sample correlation and differential peaks heatmaps
# =============================================================================
# Usage:
#   Rscript scripts/plot_heatmap.R \
#       --counts-rds results/diffbind/dba_counts.rds \
#       --sig-dir results/diffbind \
#       --output-corr-pdf results/plots/sample_correlation_heatmap.pdf \
#       --output-corr-png results/plots/sample_correlation_heatmap.png \
#       --output-diff-pdf results/plots/diff_peaks_heatmap.pdf \
#       --output-diff-png results/plots/diff_peaks_heatmap.png
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(DiffBind)
  library(pheatmap)
  library(RColorBrewer)
})

source("scripts/utils.R")

option_list <- list(
  make_option("--counts-rds", type = "character",
              help = "DBA counts object (RDS)"),
  make_option("--sig-dir", type = "character",
              help = "Directory containing .sig.tsv files"),
  make_option("--output-corr-pdf", type = "character"),
  make_option("--output-corr-png", type = "character"),
  make_option("--output-diff-pdf", type = "character"),
  make_option("--output-diff-png", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("=== Heatmap Generation ===\n")

dir.create(dirname(opt$`output-corr-pdf`), showWarnings = FALSE, recursive = TRUE)

# ── load DBA object ──────────────────────────────────────────────────────────
dba_obj <- readRDS(opt$`counts-rds`)

# ══════════════════════════════════════════════════════════════════════════════
# 1. Sample correlation heatmap
# ══════════════════════════════════════════════════════════════════════════════
cat("[1/2] Generating sample correlation heatmap...\n")

tryCatch({
  # Use DiffBind's built-in correlation plot to PDF
  pdf(opt$`output-corr-pdf`, width = 8, height = 7)
  dba.plotHeatmap(dba_obj, correlations = TRUE,
                  colScheme = "Blues",
                  main = "Sample Correlation (DRIPc-seq)")
  dev.off()

  png(opt$`output-corr-png`, width = 800, height = 700, res = 150)
  dba.plotHeatmap(dba_obj, correlations = TRUE,
                  colScheme = "Blues",
                  main = "Sample Correlation (DRIPc-seq)")
  dev.off()

  cat("  Saved:", opt$`output-corr-pdf`, "\n")
}, error = function(e) {
  cat("  WARNING: DiffBind heatmap failed:", conditionMessage(e), "\n")
  cat("  Attempting manual correlation heatmap...\n")

  tryCatch({
    # Manual approach: extract count matrix
    count_data <- dba.peakset(dba_obj, bRetrieve = TRUE)
    if (is(count_data, "GRanges")) {
      mat <- as.data.frame(mcols(count_data))
      # Select numeric columns
      mat <- mat[, sapply(mat, is.numeric), drop = FALSE]
    } else {
      mat <- as.matrix(count_data)
    }

    if (ncol(mat) > 1 && nrow(mat) > 1) {
      cor_mat <- cor(mat, use = "pairwise.complete.obs")

      pdf(opt$`output-corr-pdf`, width = 8, height = 7)
      pheatmap(cor_mat,
               color = colorRampPalette(brewer.pal(9, "Blues"))(100),
               main = "Sample Correlation (DRIPc-seq)",
               display_numbers = TRUE,
               number_format = "%.2f")
      dev.off()

      png(opt$`output-corr-png`, width = 800, height = 700, res = 150)
      pheatmap(cor_mat,
               color = colorRampPalette(brewer.pal(9, "Blues"))(100),
               main = "Sample Correlation (DRIPc-seq)",
               display_numbers = TRUE,
               number_format = "%.2f")
      dev.off()
    } else {
      placeholder_plot(opt$`output-corr-pdf`,
                       "Insufficient data for correlation heatmap")
      placeholder_plot(opt$`output-corr-png`,
                       "Insufficient data for correlation heatmap")
    }
  }, error = function(e2) {
    cat("  ERROR: Manual heatmap also failed:", conditionMessage(e2), "\n")
    placeholder_plot(opt$`output-corr-pdf`, "Correlation heatmap failed")
    placeholder_plot(opt$`output-corr-png`, "Correlation heatmap failed")
  })
})

# ══════════════════════════════════════════════════════════════════════════════
# 2. Significant differential peaks heatmap
# ══════════════════════════════════════════════════════════════════════════════
cat("[2/2] Generating significant differential peaks heatmap...\n")

tryCatch({
  # Collect all significant peaks across comparisons
  sig_files <- list.files(opt$`sig-dir`, pattern = "\\.sig\\.tsv$",
                          full.names = TRUE)

  all_sig <- data.frame()
  for (f in sig_files) {
    tmp <- read.table(f, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(tmp) > 0) {
      tmp$comparison <- gsub("\\.sig\\.tsv$", "", basename(f))
      all_sig <- rbind(all_sig, tmp)
    }
  }

  if (nrow(all_sig) == 0) {
    cat("  No significant peaks found. Creating placeholder.\n")
    placeholder_plot(opt$`output-diff-pdf`,
                     "No significant differential peaks")
    placeholder_plot(opt$`output-diff-png`,
                     "No significant differential peaks")
  } else {
    # Try to build a matrix of FC values per peak per comparison
    # Identify peak coordinates
    coord_cols <- intersect(c("seqnames", "start", "end", "Chr", "Start", "End",
                              "chrom"), colnames(all_sig))
    fc_col <- find_column(all_sig, c("Fold", "log2FoldChange", "logFC"))

    if (length(coord_cols) >= 3 && !is.null(fc_col)) {
      # Create peak ID
      chr_col <- coord_cols[1]
      st_col <- coord_cols[2]
      en_col <- coord_cols[3]
      all_sig$peak_id <- paste0(all_sig[[chr_col]], ":",
                                all_sig[[st_col]], "-",
                                all_sig[[en_col]])

      # Pivot to matrix
      fc_matrix <- tapply(all_sig[[fc_col]],
                          list(all_sig$peak_id, all_sig$comparison),
                          FUN = mean)
      fc_matrix[is.na(fc_matrix)] <- 0

      # Limit to manageable size for heatmap
      if (nrow(fc_matrix) > 500) {
        # Keep top peaks by variance
        row_var <- apply(fc_matrix, 1, var, na.rm = TRUE)
        fc_matrix <- fc_matrix[order(-row_var)[1:500], , drop = FALSE]
      }

      # Color palette
      col_palette <- colorRampPalette(
        rev(brewer.pal(11, "RdBu"))
      )(100)

      pdf(opt$`output-diff-pdf`, width = 10, height = 12)
      pheatmap(fc_matrix,
               color = col_palette,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = (nrow(fc_matrix) <= 50),
               main = "Significant Differential Peaks (log2FC)",
               fontsize = 8)
      dev.off()

      png(opt$`output-diff-png`, width = 1000, height = 1200, res = 150)
      pheatmap(fc_matrix,
               color = col_palette,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = (nrow(fc_matrix) <= 50),
               main = "Significant Differential Peaks (log2FC)",
               fontsize = 8)
      dev.off()

      cat("  Heatmap includes", nrow(fc_matrix), "peaks x",
          ncol(fc_matrix), "comparisons\n")
    } else {
      cat("  WARNING: Could not construct FC matrix. Creating placeholder.\n")
      placeholder_plot(opt$`output-diff-pdf`, "Could not build peak matrix")
      placeholder_plot(opt$`output-diff-png`, "Could not build peak matrix")
    }
  }

  cat("  Saved:", opt$`output-diff-pdf`, "\n")
  cat("  Saved:", opt$`output-diff-png`, "\n")

}, error = function(e) {
  cat("  ERROR:", conditionMessage(e), "\n")
  placeholder_plot(opt$`output-diff-pdf`, "Heatmap generation failed")
  placeholder_plot(opt$`output-diff-png`, "Heatmap generation failed")
})

cat("Heatmap generation complete.\n")
