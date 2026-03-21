#!/usr/bin/env Rscript
# =============================================================================
# run_enrichment.R — GO and KEGG enrichment for significant peak-associated genes
# =============================================================================
# Uses clusterProfiler for mouse (org.Mm.eg.db).
# If too few genes, the step is skipped and logged without error.
#
# Usage:
#   Rscript scripts/run_enrichment.R \
#       --sig-peaks results/diffbind/3day_vs_Sham.sig.tsv \
#       --orgdb org.Mm.eg.db \
#       --pval-cutoff 0.05 \
#       --min-genes 5 \
#       --output-go-bp results/enrichment/3day_vs_Sham.go_bp.tsv \
#       --output-go-mf results/enrichment/3day_vs_Sham.go_mf.tsv \
#       --output-kegg results/enrichment/3day_vs_Sham.kegg.tsv \
#       --output-go-pdf results/enrichment/3day_vs_Sham.go_bp.pdf \
#       --output-kegg-pdf results/enrichment/3day_vs_Sham.kegg.pdf \
#       --comparison "3day_vs_Sham"
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(ggplot2)
})

source("scripts/utils.R")

option_list <- list(
  make_option("--sig-peaks", type = "character",
              help = "Significant peaks TSV"),
  make_option("--orgdb", type = "character", default = "org.Mm.eg.db"),
  make_option("--pval-cutoff", type = "double", default = 0.05),
  make_option("--min-genes", type = "integer", default = 5),
  make_option("--output-go-bp", type = "character"),
  make_option("--output-go-mf", type = "character"),
  make_option("--output-kegg", type = "character"),
  make_option("--output-go-pdf", type = "character"),
  make_option("--output-kegg-pdf", type = "character"),
  make_option("--comparison", type = "character", default = "comparison")
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("=== Enrichment Analysis ===\n")
cat("Comparison:", opt$comparison, "\n")
cat("OrgDb:", opt$orgdb, "\n")

suppressPackageStartupMessages(library(opt$orgdb, character.only = TRUE))
orgdb <- get(opt$orgdb)

# Create output directories
for (f in c(opt$`output-go-bp`, opt$`output-go-mf`, opt$`output-kegg`,
            opt$`output-go-pdf`, opt$`output-kegg-pdf`)) {
  dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
}

# ── helper: write empty outputs ──────────────────────────────────────────────
write_empty <- function(reason) {
  cat("  Reason:", reason, "\n")
  for (f in c(opt$`output-go-bp`, opt$`output-go-mf`, opt$`output-kegg`)) {
    write.table(
      data.frame(Description = reason, stringsAsFactors = FALSE),
      f, sep = "\t", row.names = FALSE, quote = FALSE
    )
  }
  for (f in c(opt$`output-go-pdf`, opt$`output-kegg-pdf`)) {
    placeholder_plot(f, reason)
  }
}

# ── read significant peaks ───────────────────────────────────────────────────
sig_df <- read.table(opt$`sig-peaks`, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(sig_df) == 0) {
  cat("No significant peaks found.\n")
  write_empty("No significant peaks")
  quit(save = "no", status = 0)
}

# ── extract gene symbols ────────────────────────────────────────────────────
gene_col <- find_column(sig_df, c("gene_symbol", "SYMBOL", "Gene",
                                   "gene_name", "NEAREST_GENE"))

if (is.null(gene_col)) {
  # Try to get gene_id and convert
  id_col <- find_column(sig_df, c("gene_id", "geneId", "ENTREZID"))
  if (!is.null(id_col)) {
    entrez_ids <- unique(na.omit(sig_df[[id_col]]))
    entrez_ids <- entrez_ids[entrez_ids != "" & entrez_ids != "NA"]
  } else {
    cat("WARNING: No gene column found in sig peaks.\n")
    write_empty("No gene column in input")
    quit(save = "no", status = 0)
  }
} else {
  symbols <- unique(na.omit(sig_df[[gene_col]]))
  symbols <- symbols[symbols != "" & symbols != "NA"]

  if (length(symbols) == 0) {
    cat("No gene symbols found in significant peaks.\n")
    write_empty("No gene symbols found")
    quit(save = "no", status = 0)
  }

  # Convert gene symbols to Entrez IDs
  cat("Converting", length(symbols), "gene symbols to Entrez IDs...\n")
  id_map <- tryCatch(
    AnnotationDbi::select(orgdb,
                          keys = symbols,
                          columns = "ENTREZID",
                          keytype = "SYMBOL"),
    error = function(e) {
      cat("WARNING: Symbol to Entrez mapping failed:", conditionMessage(e), "\n")
      data.frame(SYMBOL = character(0), ENTREZID = character(0))
    }
  )
  entrez_ids <- unique(na.omit(id_map$ENTREZID))
}

cat("Entrez IDs for enrichment:", length(entrez_ids), "\n")

if (length(entrez_ids) < opt$`min-genes`) {
  msg <- paste0("Too few genes (", length(entrez_ids),
                " < min ", opt$`min-genes`, ") — skipping enrichment")
  cat(msg, "\n")
  write_empty(msg)
  quit(save = "no", status = 0)
}

# ══════════════════════════════════════════════════════════════════════════════
# GO Biological Process
# ══════════════════════════════════════════════════════════════════════════════
cat("[1/3] Running GO BP enrichment...\n")
go_bp <- tryCatch({
  enrichGO(gene = entrez_ids,
           OrgDb = orgdb,
           keyType = "ENTREZID",
           ont = "BP",
           pvalueCutoff = opt$`pval-cutoff`,
           pAdjustMethod = "BH",
           readable = TRUE)
}, error = function(e) {
  cat("  WARNING: GO BP failed:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
  go_bp_df <- as.data.frame(go_bp)
  write.table(go_bp_df, opt$`output-go-bp`, sep = "\t",
              row.names = FALSE, quote = FALSE)
  cat("  GO BP terms:", nrow(go_bp_df), "\n")

  # Plot top terms
  pdf(opt$`output-go-pdf`, width = 9, height = 7)
  n_show <- min(20, nrow(go_bp_df))
  print(dotplot(go_bp, showCategory = n_show,
                title = paste0("GO BP: ", opt$comparison)))
  dev.off()

  # Also save PNG
  png_path <- sub("\\.pdf$", ".png", opt$`output-go-pdf`)
  png(png_path, width = 900, height = 700, res = 150)
  print(dotplot(go_bp, showCategory = n_show,
                title = paste0("GO BP: ", opt$comparison)))
  dev.off()
} else {
  cat("  No significant GO BP terms.\n")
  write.table(data.frame(Description = "No significant GO BP terms"),
              opt$`output-go-bp`, sep = "\t", row.names = FALSE, quote = FALSE)
  placeholder_plot(opt$`output-go-pdf`, "No significant GO BP terms")
}

# ══════════════════════════════════════════════════════════════════════════════
# GO Molecular Function
# ══════════════════════════════════════════════════════════════════════════════
cat("[2/3] Running GO MF enrichment...\n")
go_mf <- tryCatch({
  enrichGO(gene = entrez_ids,
           OrgDb = orgdb,
           keyType = "ENTREZID",
           ont = "MF",
           pvalueCutoff = opt$`pval-cutoff`,
           pAdjustMethod = "BH",
           readable = TRUE)
}, error = function(e) {
  cat("  WARNING: GO MF failed:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(go_mf) && nrow(as.data.frame(go_mf)) > 0) {
  go_mf_df <- as.data.frame(go_mf)
  write.table(go_mf_df, opt$`output-go-mf`, sep = "\t",
              row.names = FALSE, quote = FALSE)
  cat("  GO MF terms:", nrow(go_mf_df), "\n")
} else {
  cat("  No significant GO MF terms.\n")
  write.table(data.frame(Description = "No significant GO MF terms"),
              opt$`output-go-mf`, sep = "\t", row.names = FALSE, quote = FALSE)
}

# ══════════════════════════════════════════════════════════════════════════════
# KEGG
# ══════════════════════════════════════════════════════════════════════════════
cat("[3/3] Running KEGG enrichment...\n")
kegg <- tryCatch({
  enrichKEGG(gene = entrez_ids,
             organism = "mmu",   # mouse
             pvalueCutoff = opt$`pval-cutoff`,
             pAdjustMethod = "BH")
}, error = function(e) {
  cat("  WARNING: KEGG failed:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
  kegg_df <- as.data.frame(kegg)
  write.table(kegg_df, opt$`output-kegg`, sep = "\t",
              row.names = FALSE, quote = FALSE)
  cat("  KEGG pathways:", nrow(kegg_df), "\n")

  pdf(opt$`output-kegg-pdf`, width = 9, height = 7)
  n_show <- min(20, nrow(kegg_df))
  print(dotplot(kegg, showCategory = n_show,
                title = paste0("KEGG: ", opt$comparison)))
  dev.off()

  png_path <- sub("\\.pdf$", ".png", opt$`output-kegg-pdf`)
  png(png_path, width = 900, height = 700, res = 150)
  print(dotplot(kegg, showCategory = n_show,
                title = paste0("KEGG: ", opt$comparison)))
  dev.off()
} else {
  cat("  No significant KEGG pathways.\n")
  write.table(data.frame(Description = "No significant KEGG pathways"),
              opt$`output-kegg`, sep = "\t", row.names = FALSE, quote = FALSE)
  placeholder_plot(opt$`output-kegg-pdf`, "No significant KEGG pathways")
}

cat("\nEnrichment analysis complete for:", opt$comparison, "\n")
