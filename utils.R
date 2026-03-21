#!/usr/bin/env Rscript
# =============================================================================
# utils.R — Shared utility functions for DRIPc-seq R scripts
# =============================================================================
# Source this file at the top of every R script:
#   source("scripts/utils.R")
# =============================================================================

#' Find the first matching column name from a list of candidates
#'
#' @param df  A data.frame
#' @param candidates  Character vector of possible column names
#' @return  The first matching column name, or NULL if none found
find_column <- function(df, candidates) {
  # Try exact match first
  for (col in candidates) {
    if (col %in% colnames(df)) return(col)
  }
  # Try case-insensitive match
  lower_cols <- tolower(colnames(df))
  for (col in candidates) {
    idx <- which(lower_cols == tolower(col))
    if (length(idx) > 0) return(colnames(df)[idx[1]])
  }
  return(NULL)
}


#' Create a placeholder plot with a text message
#'
#' Useful for error handling — produces a valid PDF or PNG so
#' Snakemake output requirements are satisfied.
#'
#' @param path   Output file path (PDF or PNG detected by extension)
#' @param msg    Text message to display in the plot
#' @param width  Plot width
#' @param height Plot height
placeholder_plot <- function(path, msg, width = 7, height = 5) {
  ext <- tolower(tools::file_ext(path))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  if (ext == "pdf") {
    pdf(path, width = width, height = height)
  } else if (ext == "png") {
    png(path, width = width * 100, height = height * 100, res = 100)
  } else {
    pdf(path, width = width, height = height)
  }

  plot.new()
  text(0.5, 0.5, msg, cex = 1.2, col = "grey40")
  dev.off()
}


#' Safely read a TSV file, returning an empty data.frame on error
#'
#' @param path  File path
#' @return  A data.frame (possibly empty)
safe_read_tsv <- function(path) {
  tryCatch(
    read.table(path, header = TRUE, sep = "\t",
               stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      message("WARNING: Could not read ", path, ": ", conditionMessage(e))
      data.frame()
    }
  )
}


#' Log a message with timestamp to stderr
#'
#' @param ...  Message parts (passed to paste0)
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
  message(msg)
}


#' Convert a vector of gene symbols to Entrez IDs using an OrgDb
#'
#' @param symbols   Character vector of gene symbols
#' @param orgdb     OrgDb object or its name as string
#' @return  Named character vector (names = symbols, values = Entrez IDs)
symbols_to_entrez <- function(symbols, orgdb) {
  if (is.character(orgdb)) {
    suppressPackageStartupMessages(library(orgdb, character.only = TRUE))
    orgdb <- get(orgdb)
  }
  suppressPackageStartupMessages(library(AnnotationDbi))

  id_map <- tryCatch(
    AnnotationDbi::select(orgdb,
                          keys = symbols,
                          columns = "ENTREZID",
                          keytype = "SYMBOL"),
    error = function(e) {
      message("WARNING: Symbol->Entrez mapping failed: ", conditionMessage(e))
      data.frame(SYMBOL = character(0), ENTREZID = character(0))
    }
  )

  result <- setNames(id_map$ENTREZID, id_map$SYMBOL)
  result <- result[!is.na(result)]
  return(result)
}


#' Ensure a directory exists
#'
#' @param path  Directory path
ensure_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
}
