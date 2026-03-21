#!/usr/bin/env Rscript
# =============================================================================
# deg_overlap.R — Overlap DRIPc-seq differential peaks with DEG list
# =============================================================================
# If no DEG file is provided or it does not exist, the pipeline skips this
# step without error.
#
# Usage:
#   Rscript scripts/deg_overlap.R \
#       --sig-peaks results/diffbind/3day_vs_Sham.sig.tsv \
#       --deg-file path/to/deg_table.tsv \
#       --output-summary results/overlap/3day_vs_Sham.overlap_summary.tsv \
#       --output-genes results/overlap/3day_vs_Sham.overlap_genes.tsv \
#       --output-venn results/overlap/3day_vs_Sham.venn.pdf \
#       --comparison "3day_vs_Sham"
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(VennDiagram)
})

source("scripts/utils.R")

option_list <- list(
  make_option("--sig-peaks", type = "character",
              help = "Significant peaks TSV (from DiffBind)"),
  make_option("--deg-file", type = "character",
              help = "DEG table (TSV/CSV)"),
  make_option("--output-summary", type = "character"),
  make_option("--output-genes", type = "character"),
  make_option("--output-venn", type = "character"),
  make_option("--comparison", type = "character", default = "comparison")
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("=== DEG Overlap Analysis ===\n")
cat("Comparison:", opt$comparison, "\n")

dir.create(dirname(opt$`output-summary`), showWarnings = FALSE, recursive = TRUE)

# ── read significant peaks ───────────────────────────────────────────────────
sig_df <- read.table(opt$`sig-peaks`, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, check.names = FALSE)

# Extract gene symbols from peaks
peak_gene_col <- find_column(sig_df, c("gene_symbol", "SYMBOL",
                                        "Gene", "gene_name",
                                        "NEAREST_GENE"))
if (is.null(peak_gene_col)) {
  # Try to look for any column that might have gene info
  cat("WARNING: No gene_symbol column in sig peaks. Attempting to use",
      "nearest gene annotation.\n")
  peak_genes <- character(0)
} else {
  peak_genes <- unique(na.omit(sig_df[[peak_gene_col]]))
  peak_genes <- peak_genes[peak_genes != "" & peak_genes != "NA"]
}

cat("Peak-associated genes:", length(peak_genes), "\n")

# ── read DEG file ────────────────────────────────────────────────────────────
if (!file.exists(opt$`deg-file`)) {
  cat("WARNING: DEG file not found:", opt$`deg-file`, "\n")
  cat("Skipping overlap analysis.\n")
  # Write empty outputs
  write.table(data.frame(metric = "skipped", value = "DEG file not found"),
              opt$`output-summary`, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(gene_symbol = character(0)),
              opt$`output-genes`, sep = "\t", row.names = FALSE, quote = FALSE)
  pdf(opt$`output-venn`, width = 6, height = 6)
  plot.new(); text(0.5, 0.5, "DEG file not found", cex = 1.2)
  dev.off()
  quit(save = "no", status = 0)
}

# Try both CSV and TSV
deg_df <- tryCatch(
  read.table(opt$`deg-file`, header = TRUE, sep = "\t",
             stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) {
    read.csv(opt$`deg-file`, stringsAsFactors = FALSE, check.names = FALSE)
  }
)

deg_gene_col <- find_column(deg_df, c("gene_symbol", "SYMBOL", "Gene",
                                       "gene_name", "gene_id", "GeneSymbol"))
if (is.null(deg_gene_col)) {
  cat("WARNING: Cannot find gene column in DEG file.\n")
  cat("Columns available:", paste(colnames(deg_df), collapse = ", "), "\n")
  deg_genes <- character(0)
} else {
  deg_genes <- unique(na.omit(deg_df[[deg_gene_col]]))
  deg_genes <- deg_genes[deg_genes != "" & deg_genes != "NA"]
}

cat("DEG genes:", length(deg_genes), "\n")

# ── compute overlap ──────────────────────────────────────────────────────────
overlap_genes <- intersect(peak_genes, deg_genes)
peak_only <- setdiff(peak_genes, deg_genes)
deg_only <- setdiff(deg_genes, peak_genes)

cat("Overlap genes:", length(overlap_genes), "\n")

# ── write outputs ────────────────────────────────────────────────────────────
# Summary
summary_df <- data.frame(
  metric = c("comparison", "n_peak_genes", "n_deg_genes",
             "n_overlap", "n_peak_only", "n_deg_only"),
  value = c(opt$comparison, length(peak_genes), length(deg_genes),
            length(overlap_genes), length(peak_only), length(deg_only)),
  stringsAsFactors = FALSE
)
write.table(summary_df, opt$`output-summary`, sep = "\t",
            row.names = FALSE, quote = FALSE)

# Overlap genes
if (length(overlap_genes) > 0) {
  genes_df <- data.frame(gene_symbol = overlap_genes, stringsAsFactors = FALSE)
  # Optionally merge peak and DEG info
} else {
  genes_df <- data.frame(gene_symbol = character(0))
}
write.table(genes_df, opt$`output-genes`, sep = "\t",
            row.names = FALSE, quote = FALSE)

# ── Venn diagram ─────────────────────────────────────────────────────────────
cat("Generating Venn diagram...\n")

tryCatch({
  futile.logger::flog.threshold(futile.logger::ERROR)  # suppress VennDiagram logs

  venn_list <- list(
    "DRIPc-seq peaks" = peak_genes,
    "DEGs" = deg_genes
  )

  pdf(opt$`output-venn`, width = 6, height = 6)
  venn.plot <- venn.diagram(
    x = venn_list,
    filename = NULL,
    fill = c("#E41A1C", "#377EB8"),
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 1.2,
    main = paste0("DRIPc-seq ∩ DEG: ", opt$comparison),
    main.cex = 1.3
  )
  grid::grid.draw(venn.plot)
  dev.off()

  cat("Venn diagram saved:", opt$`output-venn`, "\n")
}, error = function(e) {
  cat("WARNING: Venn diagram failed:", conditionMessage(e), "\n")
  placeholder_plot(opt$`output-venn`, "Venn diagram failed")
})

cat("DEG overlap analysis complete.\n")
