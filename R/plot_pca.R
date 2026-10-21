
#' Generate a PCA Plot for Gene Expression Data
#'
#' @description
#' Performs principal component analysis on transformed gene expression
#' data and visualizes samples using the first two principal components.
#'
#' @param expression_matrix Numeric matrix with genes as rows and
#'   samples as columns. Values should be normalized and transformed.
#' @param metadata Data frame containing sample information.
#'   Row names must match expression matrix column names.
#' @param color_by Metadata column used to color samples.
#' @param shape_by Optional metadata column used to shape samples.
#' @param top_n Number of most variable genes to use for PCA.
#' @param title_text Title displayed above the PCA plot.
#'
#' @return A ggplot object.
#'
#' @importFrom rlang .data
#' @export
plot_pca <- function(
    expression_matrix,
    metadata,
    color_by = "condition",
    shape_by = NULL,
    top_n = 500,
    title_text = "PCA of RNA-seq Samples"
) {

  # 1. Validate inputs
  if (!is.matrix(expression_matrix) ||
      !is.numeric(expression_matrix)) {
    stop("expression_matrix must be a numeric matrix.")
  }

  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame.")
  }

  if (!identical(colnames(expression_matrix),
                 rownames(metadata))) {
    stop("Sample names and order must match between data and metadata.")
  }

  if (!color_by %in% colnames(metadata)) {
    stop("color_by column not found in metadata.")
  }

  if (!is.null(shape_by) &&
      !shape_by %in% colnames(metadata)) {
    stop("shape_by column not found in metadata.")
  }

  if (!is.numeric(top_n) || length(top_n) != 1 ||
      is.na(top_n) || !is.finite(top_n) ||
      top_n < 2 || top_n != floor(top_n)) {
    stop("top_n must be a whole number greater than or equal to 2.")
  }

  if (ncol(expression_matrix) < 3) {
    stop("At least three samples are required.")
  }

  if (any(!is.finite(expression_matrix))) {
    stop("Expression matrix contains missing or non-finite values.")
  }

  # 2. Select the most variable genes
  gene_variance <- apply(expression_matrix, 1, stats::var)

  gene_variance <- sort(
    gene_variance[is.finite(gene_variance) & gene_variance > 0],
    decreasing = TRUE
  )

  if (length(gene_variance) < 2) {
    stop("At least two genes with nonzero variance are required.")
  }

  selected_genes <- names(
    utils::head(
      gene_variance,
      min(top_n, length(gene_variance))
    )
  )

  pca_matrix <- expression_matrix[selected_genes, , drop = FALSE]

  # 3. Perform PCA
  pca_result <- stats::prcomp(
    t(pca_matrix),
    center = TRUE,
    scale. = FALSE
  )

  # 4. Calculate variance explained
  percent_variance <- 100 *
    pca_result$sdev^2 / sum(pca_result$sdev^2)

  # 5. Prepare plotting data
  pca_df <- as.data.frame(pca_result$x[, 1:2, drop = FALSE])

  pca_df <- cbind(
    pca_df,
    metadata[, , drop = FALSE]
  )

  # 6. Create PCA plot
  p <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(
      x = .data$PC1,
      y = .data$PC2,
      color = .data[[color_by]]
    )
  )

  if (!is.null(shape_by)) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(shape = .data[[shape_by]]),
      size = 4
    )
  } else {
    p <- p + ggplot2::geom_point(size = 4)
  }

  p <- p +
    ggplot2::labs(
      title = title_text,
      x = paste0("PC1: ", round(percent_variance[1], 1), "% variance"),
      y = paste0("PC2: ", round(percent_variance[2], 1), "% variance"),
      color = color_by,
      shape = shape_by
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  return(p)
}
