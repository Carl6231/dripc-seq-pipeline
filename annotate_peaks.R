#!/usr/bin/env Rscript
# =============================================================================
# annotate_peaks.R — Annotate filtered peaks with overlapping transcripts
# =============================================================================
# Assigns each peak a region_class (promoter, gene_body, terminator, intergenic)
# and overlapping transcript information.
#
# Usage:
#   Rscript scripts/annotate_peaks.R \
#       --peaks results/peaks/Sham1.filtered.bed \
#       --output results/annotation/Sham1.annotated.tsv \
#       --txdb TxDb.Mmusculus.UCSC.mm10.knownGene \
#       --orgdb org.Mm.eg.db \
#       --promoter-up 2000 \
#       --promoter-down 200 \
#       --terminator-down 2000
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(ChIPseeker)
  library(AnnotationDbi)
})

source("scripts/utils.R")

# ── parse arguments ──────────────────────────────────────────────────────────
option_list <- list(
  make_option("--peaks", type = "character", help = "Input filtered BED file"),
  make_option("--output", type = "character", help = "Output annotated TSV"),
  make_option("--txdb", type = "character",
              default = "TxDb.Mmusculus.UCSC.mm10.knownGene"),
  make_option("--orgdb", type = "character", default = "org.Mm.eg.db"),
  make_option("--promoter-up", type = "integer", default = 2000),
  make_option("--promoter-down", type = "integer", default = 200),
  make_option("--terminator-down", type = "integer", default = 2000)
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("=== Peak Annotation ===\n")
cat("Input peaks:", opt$peaks, "\n")
cat("TxDb:", opt$txdb, "\n")
cat("OrgDb:", opt$orgdb, "\n")

# ── load annotation databases ────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(opt$txdb, character.only = TRUE)
  library(opt$orgdb, character.only = TRUE)
})

txdb <- get(opt$txdb)
orgdb <- get(opt$orgdb)

# ── read peaks ───────────────────────────────────────────────────────────────
peaks_df <- read.table(opt$peaks, header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE, fill = TRUE)

if (nrow(peaks_df) == 0) {
  cat("WARNING: No peaks in input file. Writing empty output.\n")
  empty_df <- data.frame(
    chrom = character(0), start = integer(0), end = integer(0),
    width = integer(0), strand = character(0),
    transcript_id = character(0), gene_id = character(0),
    gene_symbol = character(0), transcript_biotype = character(0),
    region_class = character(0), stringsAsFactors = FALSE
  )
  dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)
  write.table(empty_df, opt$output, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Empty annotation written to:", opt$output, "\n")
  quit(save = "no", status = 0)
}

# Ensure at least 3 columns
colnames(peaks_df)[1:min(ncol(peaks_df), 6)] <-
  c("chrom", "start", "end", "name", "score", "strand")[1:min(ncol(peaks_df), 6)]

peaks_gr <- GRanges(
  seqnames = peaks_df$chrom,
  ranges = IRanges(start = as.integer(peaks_df$start) + 1,  # BED is 0-based
                   end = as.integer(peaks_df$end)),
  strand = if ("strand" %in% colnames(peaks_df) && any(peaks_df$strand %in% c("+", "-")))
             peaks_df$strand else "*"
)

cat("Loaded", length(peaks_gr), "peaks\n")

# ── annotate with ChIPseeker ────────────────────────────────────────────────
cat("Annotating peaks with ChIPseeker...\n")
anno <- annotatePeak(
  peaks_gr, TxDb = txdb,
  tssRegion = c(-opt$`promoter-up`, opt$`promoter-down`),
  level = "transcript"
)

anno_df <- as.data.frame(anno)

# ── build region_class ───────────────────────────────────────────────────────
# ChIPseeker annotation field is in 'annotation'
assign_region_class <- function(chipseeker_annotation) {
  ann <- tolower(chipseeker_annotation)
  ifelse(grepl("promoter", ann), "promoter",
  ifelse(grepl("exon|intron|utr|cds|coding", ann), "gene_body",
  ifelse(grepl("downstream|3'", ann), "terminator",
  ifelse(grepl("intergenic|distal", ann), "intergenic",
         "gene_body"))))
}

anno_df$region_class <- assign_region_class(anno_df$annotation)

# ── get gene symbols from orgdb ──────────────────────────────────────────────
cat("Mapping gene IDs to symbols...\n")
if ("geneId" %in% colnames(anno_df)) {
  gene_ids <- unique(na.omit(anno_df$geneId))
  if (length(gene_ids) > 0) {
    sym_map <- tryCatch({
      AnnotationDbi::select(orgdb,
                            keys = gene_ids,
                            columns = c("SYMBOL"),
                            keytype = "ENTREZID")
    }, error = function(e) {
      cat("  WARNING: Could not map gene symbols:", conditionMessage(e), "\n")
      data.frame(ENTREZID = character(0), SYMBOL = character(0))
    })
    anno_df <- merge(anno_df, sym_map,
                     by.x = "geneId", by.y = "ENTREZID",
                     all.x = TRUE, sort = FALSE)
    anno_df$gene_symbol <- anno_df$SYMBOL
  } else {
    anno_df$gene_symbol <- NA
  }
} else {
  anno_df$gene_symbol <- NA
}

# ── assemble output columns ─────────────────────────────────────────────────
# Map column names to output schema
out_df <- data.frame(
  chrom       = anno_df$seqnames,
  start       = anno_df$start - 1,  # back to 0-based BED
  end         = anno_df$end,
  width       = anno_df$width,
  strand      = as.character(anno_df$strand),
  transcript_id    = ifelse(is.null(anno_df$transcriptId),
                            NA, anno_df$transcriptId),
  gene_id          = ifelse(is.null(anno_df$geneId),
                            NA, anno_df$geneId),
  gene_symbol      = anno_df$gene_symbol,
  transcript_biotype = ifelse("transcriptBiotype" %in% colnames(anno_df),
                              anno_df$transcriptBiotype, NA),
  region_class     = anno_df$region_class,
  chipseeker_annotation = anno_df$annotation,
  stringsAsFactors = FALSE
)

# ── write output ─────────────────────────────────────────────────────────────
dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)
write.table(out_df, opt$output, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Annotation complete:", nrow(out_df), "peaks annotated\n")
cat("Region class summary:\n")
print(table(out_df$region_class))
cat("Output:", opt$output, "\n")
