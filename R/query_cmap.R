#' Query CMap/LINCS database with a gene signature
#'
#' @param signature A list with upgenes and downgenes (from build_signature)
#' @param db Database to query, either "lincs" or "lincs2" (default "lincs")
#' @param n_top Number of top results to return (default 100)
#' @param tau Logical, whether to compute tau score (default FALSE)
#' @param workers Number of parallel workers for database query (default 1)
#' @param add_annotations Logical, whether to merge LINCS compound annotations (default TRUE)
#'
#' @return A dataframe of top drug candidates with scores and LINCS annotations
#' @export
query_cmap <- function(signature, db = "lincs", n_top = 100, tau = FALSE, workers = 1, add_annotations = TRUE) #' Query CMap/LINCS database with a gene signature
  #'
  #' @param signature A list with upgenes and downgenes (from build_signature)
  #' @param db Database to query, either "lincs" or "lincs2" (default "lincs")
  #' @param n_top Number of top results to return (default 100)
  #' @param tau Logical, whether to compute tau score (default FALSE)
  #' @param workers Number of parallel workers for database query (default 1)
  #' @param add_annotations Logical, whether to merge LINCS compound annotations (default TRUE)
  #'
  #' @return A dataframe of top drug candidates with scores and LINCS annotations
  #' @export
  query_cmap <- function(signature, db = "lincs", n_top = 100, tau = FALSE, workers = 1, add_annotations = TRUE) {

    # Check signature structure
    if (!is.list(signature) || !all(c("upgenes", "downgenes") %in% names(signature))) {
      stop("`signature` must be a list containing `upgenes` and `downgenes`.")
    }

    up <- as.character(signature$upgenes)
    down <- as.character(signature$downgenes)

    # Convert gene symbols to Entrez IDs if necessary
    .convert_genes <- function(genes) {
      if (length(genes) == 0) return(character(0))
      genes <- as.character(genes)
      genes <- unique(genes[!is.na(genes) & genes != ""])
      if (length(genes) == 0) return(character(0))
      is_numeric <- grepl("^[0-9]+$", genes)
      result <- genes
      symbols_to_convert <- genes[!is_numeric]
      if (length(symbols_to_convert) > 0) {
        converted <- as.character(AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db,
          keys = symbols_to_convert,
          column = "ENTREZID",
          keytype = "SYMBOL",
          multiVals = "first"
        ))
        result[!is_numeric] <- converted
      }
      result <- unique(result[!is.na(result) & result != ""])
      return(result)
    }

    up <- .convert_genes(up)
    down <- .convert_genes(down)

    if (length(up) == 0) {
      warning("No valid up-regulated genes after conversion.")
    }
    if (length(down) == 0) {
      warning("No valid down-regulated genes after conversion.")
    }
    if (length(up) == 0 && length(down) == 0) {
      stop("No valid genes after Entrez ID conversion. Please check gene symbols or Entrez IDs.")
    }

    # Build query signature
    qsig <- signatureSearch::qSig(
      query = list(upset = up, downset = down),
      gess_method = "LINCS",
      refdb = db
    )

    # Run LINCS query
    result <- tryCatch(
      {
        signatureSearch::gess_lincs(
          qSig = qsig,
          sortby = "NCS",
          tau = tau,
          addAnnotations = FALSE,
          workers = workers
        )
      },
      error = function(e) {
        stop(
          "CMap/LINCS query failed. Please check whether the LINCS reference database is installed and accessible. Original error: ",
          e$message
        )
      }
    )

    result_df <- signatureSearch::result(result)
    result_df <- head(result_df, n_top)
    result_df <- as.data.frame(result_df)

    # Preserve original ranking
    result_df$.rank <- seq_len(nrow(result_df))

    # Add LINCS compound annotations
    if (isTRUE(add_annotations)) {

      if (!requireNamespace("signatureSearchData", quietly = TRUE)) {
        warning("Package `signatureSearchData` is not installed. Returning unannotated CMap results.")
        return(result_df)
      }

      result_df <- tryCatch(
        {
          annot_env <- new.env()
          data("lincs_pert_info", package = "signatureSearchData", envir = annot_env)

          if (!exists("lincs_pert_info", envir = annot_env)) {
            stop("`lincs_pert_info` could not be loaded from signatureSearchData.")
          }

          lincs_pert_info <- as.data.frame(annot_env$lincs_pert_info)

          required_cols <- c("pert_id", "pert_iname")
          missing_cols <- setdiff(required_cols, colnames(lincs_pert_info))
          if (length(missing_cols) > 0) {
            stop("Required columns missing from `lincs_pert_info`: ",
                 paste(missing_cols, collapse = ", "))
          }

          # Match by both pert_iname and pert_id
          annot_by_name <- lincs_pert_info
          annot_by_name$join_key <- annot_by_name$pert_iname

          annot_by_id <- lincs_pert_info
          annot_by_id$join_key <- annot_by_id$pert_id

          annot_long <- rbind(annot_by_name, annot_by_id)
          annot_long <- annot_long[!is.na(annot_long$join_key) & annot_long$join_key != "", , drop = FALSE]
          annot_long <- annot_long[!duplicated(annot_long$join_key), , drop = FALSE]

          # Remove overlapping columns
          overlapping_cols <- intersect(
            colnames(result_df),
            setdiff(colnames(annot_long), "join_key")
          )
          if (length(overlapping_cols) > 0) {
            annot_long <- annot_long[, !colnames(annot_long) %in% overlapping_cols, drop = FALSE]
          }

          merged_df <- merge(
            result_df,
            annot_long,
            by.x = "pert",
            by.y = "join_key",
            all.x = TRUE,
            sort = FALSE
          )

          # Restore original ranking
          merged_df <- merged_df[order(merged_df$.rank), ]
          merged_df$.rank <- NULL
          merged_df
        },
        error = function(e) {
          warning("LINCS annotation failed. Returning unannotated CMap results. Original error: ", e$message)
          result_df$.rank <- NULL
          result_df
        }
      )
    } else {
      result_df$.rank <- NULL
    }

    return(result_df)
  }
