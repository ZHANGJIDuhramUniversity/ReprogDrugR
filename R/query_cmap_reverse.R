#' Query CMap/LINCS database in reverse to identify harmful compounds
#'
#' @param signature A list with upgenes and downgenes (from build_signature)
#' @param db Database to query, either "lincs" or "lincs2" (default "lincs")
#' @param n_top Number of top results to return (default 100)
#'
#' @return A dataframe of compounds that would worsen beta cell state
#' @export
query_cmap_reverse <- function(signature,
                               db = "lincs",
                               n_top = 100) {

  # 反转签名：升调变降调，降调变升调
  reversed_signature <- list(
    upgenes = signature$downgenes,
    downgenes = signature$upgenes
  )

  # 构建查询
  qsig <- signatureSearch::qSig(
    query = list(
      upset = reversed_signature$upgenes,
      downset = reversed_signature$downgenes
    ),
    gess_method = "LINCS",
    refdb = db
  )

  # 查询数据库
  result <- signatureSearch::gess_lincs(
    qSig = qsig,
    sortby = "NCS",
    tau = TRUE
  )

  # 提取结果
  result_df <- signatureSearch::result(result)

  return(head(result_df, n_top))
}
