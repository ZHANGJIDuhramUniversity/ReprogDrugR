#' Build a gene signature from DEG results
#'
#' @param deg_file Path to the DEG CSV file
#' @param lfc_col Column name for log2FoldChange
#' @param fdr_col Column name for FDR/padj
#' @param lfc_cutoff Minimum absolute log2FoldChange (default 0.5)
#' @param fdr_cutoff Maximum FDR (default 0.05)
#'
#' @return A list with upgenes and downgenes
#' @export
build_signature <- function(deg_file,
                            lfc_col = "log2FoldChange",
                            fdr_col = "padj",
                            lfc_cutoff = 0.5,
                            fdr_cutoff = 0.05) {

  # 读取文件
  deg <- read.csv(deg_file)

  # 筛选显著基因
  deg <- deg[!is.na(deg[[fdr_col]]), ]
  deg <- deg[deg[[fdr_col]] < fdr_cutoff, ]
  deg <- deg[abs(deg[[lfc_col]]) > lfc_cutoff, ]

  # 拆分上调和下调
  upgenes <- deg$gene[deg[[lfc_col]] > 0]
  downgenes <- deg$gene[deg[[lfc_col]] < 0]

  return(list(
    upgenes = upgenes,
    downgenes = downgenes
  ))
}
