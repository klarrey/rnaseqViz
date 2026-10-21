#' Generate a Volcano Plot for Differential Gene Expression
#'
#' Creates a volcano plot from differential expression results,
#' allowing independent log2 fold-change thresholds for
#' upregulated and downregulated genes.
#'
#' @param df A data frame containing differential expression results.
#' @param title_text Character string specifying the plot title.
#' @param up_thr Positive log2 fold-change threshold.
#' @param down_thr Negative log2 fold-change threshold.
#' @param p_thr Statistical significance threshold for raw p-values.
#' @param top_n Number of top genes to label in each direction.
#'
#' @return A ggplot object.
#'
#' @export
#' @importFrom rlang .data


plot_volcano <- function(df, title_text,
                         up_thr = 1, down_thr = -1, p_thr = 0.05,
                         top_n = 10) {

  # Validate input data
  if (!is.data.frame(df)) {
    stop("df must be a data frame.", call. = FALSE)
  }

  # Check required columns
  required_cols <- c(
    "log2FoldChange",
    "pvalue",
    "gene_name"
  )

  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    stop(
      paste(
        "Missing required columns:",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Validate thresholds
  if (!is.numeric(up_thr) || length(up_thr) != 1 ||
      is.na(up_thr) || !is.finite(up_thr) || up_thr <= 0) {
    stop("up_thr must be a positive number.", call. = FALSE)
  }

  if (!is.numeric(down_thr) || length(down_thr) != 1 ||
      is.na(down_thr) || !is.finite(down_thr) || down_thr >= 0) {
    stop("down_thr must be a negative number.", call. = FALSE)
  }

  if (!is.numeric(p_thr) || length(p_thr) != 1 ||
      is.na(p_thr) || !is.finite(p_thr) ||
      p_thr <= 0 || p_thr >= 1) {
    stop("p_thr must be between 0 and 1.", call. = FALSE)
  }

  if (!is.numeric(top_n) || length(top_n) != 1 ||
      is.na(top_n) || !is.finite(top_n) ||
      top_n < 0 || top_n != floor(top_n)) {
    stop("top_n must be a non-negative whole number.",
         call. = FALSE)
  }

  # Validate thresholds
  stopifnot(
    up_thr > 0,
    down_thr < 0,
    p_thr > 0,
    p_thr < 1
  )

  # Prepare differential expression results
  df <- dplyr::mutate(
    df,
    log2FoldChange = as.numeric(.data$log2FoldChange),
    pvalue = as.numeric(.data$pvalue),
    negLogP = -log10(.data$pvalue),

    status = dplyr::case_when(
      .data$log2FoldChange >= up_thr &
        .data$pvalue < p_thr ~ "Up",

      .data$log2FoldChange <= down_thr &
        .data$pvalue < p_thr ~ "Down",

      TRUE ~ "NS"
    )
  )

  # Count genes in each category
  total_n <- nrow(df)

  up_n <- sum(df$status == "Up", na.rm = TRUE)

  down_n <- sum(df$status == "Down", na.rm = TRUE)


  #Exclude genes without valid symbols from labeling
  # Select top 10 upregulated genes
  top_up <- df |>
    dplyr::filter(
      .data$status == "Up",
      !is.na(.data$gene_name),
      .data$gene_name != ""
    ) |>
    dplyr::arrange(
      .data$pvalue,
      dplyr::desc(.data$log2FoldChange)
    ) |>
    dplyr::slice_head(n = top_n)


  # Select top 10 downregulated genes
  top_down <- df |>
    dplyr::filter(
      .data$status == "Down",
      !is.na(.data$gene_name),
      .data$gene_name != ""
    ) |>
    dplyr::arrange(
      .data$pvalue,
      .data$log2FoldChange
    ) |>
    dplyr::slice_head(n = top_n)

  # Combine genes for labeling
  label_df <- dplyr::bind_rows(
    top_up,
    top_down
  )

  # Determine symmetric x-axis limits
  max_abs <- max(
    abs(df$log2FoldChange),
    na.rm = TRUE
  )

  # Generate volcano plot
  g <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$log2FoldChange, y = .data$negLogP)
  ) +

    ggplot2::geom_point(
      ggplot2::aes(color = .data$status),
      alpha = 0.7,
      size = 1.8
    ) +

    ggplot2::scale_color_manual(
      values = c(
        "Up" = "red",
        "Down" = "blue",
        "NS" = "grey70"
      ),
      breaks = c("Up", "Down", "NS")
    ) +

    ggplot2::geom_vline(
      xintercept = c(down_thr, up_thr),
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.5
    ) +

    ggplot2::geom_hline(
      yintercept = -log10(p_thr),
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.5
    ) +

    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = .data$gene_name),
      size = 3.5,
      max.overlaps = Inf,
      box.padding = 0.8,
      point.padding = 0.5,
      seed = 123
    ) +

    ggplot2::labs(
      title = title_text,
      subtitle = paste0(
        "Total: ", total_n,
        " | Up: ", up_n,
        " | Down: ", down_n,
        " | p < ", p_thr
      ),
      x = "log2 fold change",
      y = expression(-log[10](P.Value)),
      color = ""
    ) +

    ggplot2::coord_cartesian(
      xlim = c(-max_abs, max_abs)
    ) +

    ggplot2::theme_bw(base_size = 15) +

    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "right"
    )

  return(g)
}
