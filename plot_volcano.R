#!/usr/bin/env Rscript
# =============================================================================
# plot_volcano.R — Volcano plot for differential peak analysis results
# =============================================================================
# Usage:
#   Rscript scripts/plot_volcano.R \
#       --input results/diffbind/3day_vs_Sham.all.tsv \
#       --output-pdf results/plots/3day_vs_Sham.volcano.pdf \
#       --output-png results/plots/3day_vs_Sham.volcano.png \
#       --fc-threshold 1.5 \
#       --pval-threshold 0.05 \
#       --title "3day_vs_Sham"
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
})

source("scripts/utils.R")

option_list <- list(
  make_option("--input", type = "character", help = "Input all-peaks TSV from DiffBind"),
  make_option("--output-pdf", type = "character", help = "Output PDF path"),
  make_option("--output-png", type = "character", help = "Output PNG path"),
  make_option("--fc-threshold", type = "double", default = 1.5),
  make_option("--pval-threshold", type = "double", default = 0.05),
  make_option("--title", type = "character", default = "Volcano Plot")
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("=== Volcano Plot ===\n")
cat("Input:", opt$input, "\n")

# ── read data ────────────────────────────────────────────────────────────────
df <- read.table(opt$input, header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(df) == 0) {
  cat("WARNING: Empty input. Generating placeholder plot.\n")
  pdf(opt$`output-pdf`, width = 7, height = 6)
  plot.new()
  text(0.5, 0.5, paste("No data for", opt$title), cex = 1.5)
  dev.off()
  png(opt$`output-png`, width = 700, height = 600)
  plot.new()
  text(0.5, 0.5, paste("No data for", opt$title), cex = 1.5)
  dev.off()
  quit(save = "no", status = 0)
}

# ── identify fold-change and p-value columns ─────────────────────────────────
fc_col <- find_column(df, c("Fold", "log2FoldChange", "logFC", "Log2FC"))
pval_col <- find_column(df, c("p.value", "p-value", "pvalue", "Pval", "PValue"))
fdr_col <- find_column(df, c("FDR", "padj", "q.value", "qvalue"))

if (is.null(fc_col) || is.null(pval_col)) {
  cat("ERROR: Could not find fold-change or p-value columns.\n")
  cat("Available columns:", paste(colnames(df), collapse = ", "), "\n")
  # Still create empty outputs
  pdf(opt$`output-pdf`, width = 7, height = 6)
  plot.new()
  text(0.5, 0.5, "Could not identify FC/pval columns", cex = 1.2)
  dev.off()
  png(opt$`output-png`, width = 700, height = 600)
  plot.new()
  text(0.5, 0.5, "Could not identify FC/pval columns", cex = 1.2)
  dev.off()
  quit(save = "no", status = 0)
}

log2fc_thresh <- log2(opt$fc_threshold)

# ── classify peaks ───────────────────────────────────────────────────────────
df$neg_log10p <- -log10(df[[pval_col]])
df$significance <- "NS"
df$significance[df[[fc_col]] >= log2fc_thresh & df[[pval_col]] < opt$pval_threshold] <- "Up"
df$significance[df[[fc_col]] <= -log2fc_thresh & df[[pval_col]] < opt$pval_threshold] <- "Down"

n_up <- sum(df$significance == "Up")
n_down <- sum(df$significance == "Down")
n_ns <- sum(df$significance == "NS")

cat("Up:", n_up, " Down:", n_down, " NS:", n_ns, "\n")

# ── plot ─────────────────────────────────────────────────────────────────────
colors <- c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey60")

p <- ggplot(df, aes(x = .data[[fc_col]], y = neg_log10p, color = significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = colors,
                     labels = c(paste0("Down (", n_down, ")"),
                                paste0("NS (", n_ns, ")"),
                                paste0("Up (", n_up, ")"))) +
  geom_vline(xintercept = c(-log2fc_thresh, log2fc_thresh),
             linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(opt$pval_threshold),
             linetype = "dashed", color = "grey40") +
  labs(
    title = gsub("_", " ", opt$title),
    subtitle = paste0("FC threshold: ", opt$fc_threshold,
                      " | p < ", opt$pval_threshold),
    x = "log2(Fold Change)",
    y = "-log10(p-value)",
    color = "Significance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

# ── save ─────────────────────────────────────────────────────────────────────
dir.create(dirname(opt$`output-pdf`), showWarnings = FALSE, recursive = TRUE)

ggsave(opt$`output-pdf`, p, width = 7, height = 6)
ggsave(opt$`output-png`, p, width = 7, height = 6, dpi = 300)

cat("Saved:", opt$`output-pdf`, "\n")
cat("Saved:", opt$`output-png`, "\n")
