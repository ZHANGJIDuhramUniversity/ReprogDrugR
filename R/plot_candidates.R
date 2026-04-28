#' Plot drug candidates for beta cell reprogramming
#'
#' @param forward_results Dataframe from score_reprogramming()
#' @param reverse_results Dataframe from flag_reprogramming_blockers()
#' @param top_n Number of top candidates to label (default 20)
#' @param type Type of plot: "bubble", "bar", "scatter", or "heatmap"
#'
#' @return A ggplot object
#' @export
plot_candidates <- function(forward_results,
                            reverse_results = NULL,
                            top_n = 20,
                            type = "bubble") {

  df <- head(forward_results, top_n)

  pert <- reprogramming_score <- score_ncs <- score_pathway <- flag <- NULL

  if (type == "bubble") {
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = score_ncs,
      y = reprogramming_score,
      size = score_pathway,
      color = flag,
      label = pert
    )) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_text(size = 3, vjust = -1) +
      ggplot2::scale_color_manual(values = c(
        "safe" = "#2ecc71",
        "blocker" = "#e67e22",
        "toxic" = "#e74c3c"
      )) +
      ggplot2::labs(
        title = "Beta Cell Reprogramming Drug Candidates",
        x = "CMap NCS Score",
        y = "Reprogramming Score",
        size = "Pathway Score",
        color = "Safety"
      ) +
      ggplot2::theme_bw()

  } else if (type == "bar") {
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = reorder(pert, reprogramming_score),
      y = reprogramming_score,
      fill = flag
    )) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values = c(
        "safe" = "#2ecc71",
        "blocker" = "#e67e22",
        "toxic" = "#e74c3c"
      )) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Top Drug Candidates",
        x = "Drug",
        y = "Reprogramming Score",
        fill = "Safety"
      ) +
      ggplot2::theme_bw()

  } else if (type == "scatter") {
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = score_ncs,
      y = reprogramming_score,
      color = flag,
      label = pert
    )) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_text(size = 3, vjust = -1) +
      ggplot2::scale_color_manual(values = c(
        "safe" = "#2ecc71",
        "blocker" = "#e67e22",
        "toxic" = "#e74c3c"
      )) +
      ggplot2::labs(
        title = "Forward vs Reprogramming Score",
        x = "CMap NCS Score",
        y = "Reprogramming Score",
        color = "Safety"
      ) +
      ggplot2::theme_bw()
  }

  return(p)
}
