#' Score drug candidates for beta cell reprogramming potential
#'
#' @param forward_results Dataframe from query_cmap()
#' @param pathway_genes Named list of beta cell reprogramming pathway genes
#' @param weight_ncs Weight for CMap NCS score (default 0.3)
#' @param weight_pathway Weight for pathway score (default 0.3)
#' @param weight_literature Weight for literature score (default 0.4)
#'
#' @return A dataframe of drug candidates with comprehensive reprogramming scores
#' @export
score_reprogramming <- function(forward_results,
                                pathway_genes,
                                weight_ncs = 0.3,
                                weight_pathway = 0.3,
                                weight_literature = 0.4) {

  results <- forward_results

  # 1. 标准化NCS分数到0-1
  results$score_ncs <- (results$NCS - min(results$NCS)) /
    (max(results$NCS) - min(results$NCS))

  # 2. 计算通路分数
  # 检查每个药物是否与关键通路基因有关联
  pathway_gene_list <- unlist(pathway_genes)

  results$score_pathway <- sapply(results$pert, function(drug) {
    # 在通路基因里搜索药物名（简单匹配）
    matched <- any(grepl(drug, pathway_gene_list, ignore.case = TRUE))
    return(as.numeric(matched))
  })

  # 3. 文献分数（默认为0，用户可以手动更新）
  results$score_literature <- 0

  # 4. 计算综合分数
  results$reprogramming_score <-
    weight_ncs * results$score_ncs +
    weight_pathway * results$score_pathway +
    weight_literature * results$score_literature

  # 5. 按综合分数排序
  results <- results[order(results$reprogramming_score, decreasing = TRUE), ]

  return(results)
}
