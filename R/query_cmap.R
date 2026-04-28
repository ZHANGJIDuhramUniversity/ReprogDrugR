#' Query CMap/LINCS database with a gene signature
#'
#' @param signature A list with upgenes and downgenes (from build_signature)
#' @param db Database to query, either "lincs" or "lincs2" (default "lincs")
#' @param n_top Number of top results to return (default 100)
#'
#' @return A dataframe of top drug candidates with scores
#' @export
query_cmap <- function(signature,
                       db = "lincs",
                       n_top = 100) {

  # 构建查询
  qsig <- signatureSearch::qSig(
    query = list(
      upset = signature$upgenes,
      downset = signature$downgenes
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

  # 返回前n个
  return(head(result_df, n_top))
}
