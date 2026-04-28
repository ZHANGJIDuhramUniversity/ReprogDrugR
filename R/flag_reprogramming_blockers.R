#' Flag compounds that block beta cell reprogramming
#'
#' @param forward_results Dataframe from query_cmap()
#' @param reverse_results Dataframe from query_cmap_reverse()
#' @param toxic_drugs Character vector of known beta cell toxic drugs
#'
#' @return A dataframe with a 'flag' column marking dangerous compounds
#' @export
flag_reprogramming_blockers <- function(forward_results,
                                        reverse_results,
                                        toxic_drugs = c("streptozotocin",
                                                        "alloxan",
                                                        "STZ")) {

  results <- forward_results

  # 初始化flag列
  results$flag <- "safe"
  results$flag_reason <- ""

  # 1. 标记反向查询中也出现的药物（矛盾药物）
  reverse_drugs <- reverse_results$pert
  contradictory <- results$pert %in% reverse_drugs
  results$flag[contradictory] <- "blocker"
  results$flag_reason[contradictory] <- "appears in reverse query"

  # 2. 标记毒性黑名单药物
  toxic <- sapply(results$pert, function(drug) {
    any(grepl(drug, toxic_drugs, ignore.case = TRUE))
  })
  results$flag[toxic] <- "toxic"
  results$flag_reason[toxic] <- "known beta cell toxin"

  # 3. 安全药物排在前面
  results <- results[order(results$flag == "safe", decreasing = TRUE), ]

  return(results)
}
