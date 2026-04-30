#' Build a gene signature from differential expression results
#'
#' Reads a DEG results file and extracts up- and down-regulated gene
#' signatures for input to query_cmap(). Supports flexible filtering
#' by log2FoldChange and FDR thresholds, with independent control of
#' up/down gene counts via n_up and n_down parameters.
#'
#' @param deg_file Path to the DEG CSV file. Must contain columns for
#'   gene symbols, log2FoldChange, and FDR/padj.
#' @param gene_col Column name for gene symbols (default "gene_symbol").
#' @param lfc_col Column name for log2FoldChange (default "log2FoldChange").
#' @param fdr_col Column name for adjusted p-value (default "padj").
#' @param lfc_cutoff Minimum absolute log2FoldChange threshold
#'   (default 1, i.e. 2-fold change).
#' @param fdr_cutoff Maximum FDR threshold (default 0.05).
#' @param top_n Maximum number of genes per direction when n_up and
#'   n_down are not specified (default 150).
#' @param n_up Maximum number of up-regulated genes to include
#'   (default NULL, uses top_n). Use to reduce query size for speed.
#' @param n_down Maximum number of down-regulated genes to include
#'   (default NULL, uses top_n). Use to reduce query size for speed.
#'
#' @return A named list compatible with query_cmap():
#'   \itemize{
#'     \item \code{upgenes}: Top up-regulated gene symbols
#'     \item \code{downgenes}: Top down-regulated gene symbols
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' # Standard run
#' sig <- build_signature(
#'   deg_file = "/home/rstudio/data/results/deg_T2D_vs_ND_pseudobulk.csv"
#' )
#'
#' # Smaller query for faster runtime
#' sig_small <- build_signature(
#'   deg_file = "/home/rstudio/data/results/deg_T2D_vs_ND_pseudobulk.csv",
#'   n_up     = 50,
#'   n_down   = 50
#' )
#'
#' cat("Up genes:", length(sig$upgenes), "\n")
#' cat("Down genes:", length(sig$downgenes), "\n")
#' }
build_signature <- function(deg_file,
                            gene_col   = "gene_symbol",
                            lfc_col    = "log2FoldChange",
                            fdr_col    = "padj",
                            lfc_cutoff = 1,
                            fdr_cutoff = 0.05,
                            top_n      = 150,
                            n_up       = NULL,
                            n_down     = NULL) {

  # ---------- Input validation ----------
  if (!file.exists(deg_file))
    stop("deg_file not found: ", deg_file)

  # ---------- Read DEG file ----------
  deg <- read.csv(deg_file, stringsAsFactors = FALSE)

  # Check required columns
  for (col in c(gene_col, lfc_col, fdr_col)) {
    if (!col %in% colnames(deg))
      stop(sprintf("Column '%s' not found. Available: %s",
                   col, paste(colnames(deg), collapse = ", ")))
  }

  # ---------- Filter: remove NA gene symbols ----------
  deg <- deg[!is.na(deg[[gene_col]]) & deg[[gene_col]] != "", ]

  # ---------- Filter: remove NA FDR ----------
  deg <- deg[!is.na(deg[[fdr_col]]), ]

  # ---------- Filter: FDR and LFC thresholds ----------
  deg <- deg[deg[[fdr_col]] < fdr_cutoff, ]
  deg <- deg[abs(deg[[lfc_col]]) > lfc_cutoff, ]

  if (nrow(deg) == 0)
    stop(sprintf(
      "No genes passed filters (FDR < %.2f, |LFC| > %.1f). ",
      fdr_cutoff, lfc_cutoff,
      "Try relaxing lfc_cutoff or fdr_cutoff."
    ))

  # ---------- Resolve n_up and n_down ----------
  # n_up / n_down override top_n if specified
  n_up_final   <- if (!is.null(n_up))   n_up   else top_n
  n_down_final <- if (!is.null(n_down)) n_down else top_n

  # ---------- Extract up-regulated genes ----------
  up       <- deg[deg[[lfc_col]] > 0, ]
  up       <- up[order(up[[lfc_col]], decreasing = TRUE), ]
  upgenes  <- head(up[[gene_col]], n_up_final)

  # ---------- Extract down-regulated genes ----------
  down      <- deg[deg[[lfc_col]] < 0, ]
  down      <- down[order(down[[lfc_col]], decreasing = FALSE), ]
  downgenes <- head(down[[gene_col]], n_down_final)

  # ---------- Warn if very few genes ----------
  if (length(upgenes) < 10)
    warning(sprintf(
      "Only %d up-regulated genes found. Consider relaxing lfc_cutoff.",
      length(upgenes)
    ))
  if (length(downgenes) < 10)
    warning(sprintf(
      "Only %d down-regulated genes found. Consider relaxing lfc_cutoff.",
      length(downgenes)
    ))

  # ---------- Report ----------
  message(sprintf(
    "Signature built | Up: %d genes | Down: %d genes | LFC > %.1f | FDR < %.2f",
    length(upgenes), length(downgenes), lfc_cutoff, fdr_cutoff
  ))

  return(list(
    upgenes   = upgenes,
    downgenes = downgenes
  ))
}
