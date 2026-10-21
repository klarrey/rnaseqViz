
#' Generate a Gene Expression Heatmap
#'
#' @description
#' Creates a heatmap from transformed gene expression data.
#' Genes can be selected from differential expression results
#' or provided as a custom list.
#'
#' @param expression_matrix Numeric matrix with genes as rows
#'   and samples as columns. Values should be transformed.
#' @param metadata Data frame with samples as rows.
#'   Row names must match expression matrix column names.
#' @param de_results Optional differential expression results
#'   data frame containing gene IDs, log2FoldChange and pvalue.
#' @param gene_id_col Column containing gene IDs in de_results.
#' @param gene_name_col Column containing gene symbols or
#'   gene names in de_results, used to label heatmap rows.
#'   Default is "gene_name".
#' @param gene_list Optional vector of gene IDs to display.
#' @param top_n Number of top DE genes to select.
#' @param p_thr Raw p-value significance threshold.
#' @param up_thr Positive log2 fold-change threshold.
#' @param down_thr Negative log2 fold-change threshold.
#' @param scale_rows Apply row-wise Z-score normalization.
#' @param cluster_rows Cluster genes.
#' @param cluster_cols Cluster samples.
#' @param annotation_cols Metadata columns used for annotation.
#' @param show_rownames Display gene identifiers.
#' @param title_text Heatmap title.
#'
#' @return A pheatmap object, returned invisibly.
#'
#' @export
plot_heatmap <- function(
    expression_matrix,
    metadata,
    de_results = NULL,
    gene_id_col = "entrez_id",
    gene_name_col = "gene_name",
    gene_list = NULL,
    top_n = 10,
    p_thr = 0.05,
    up_thr = 1,
    down_thr = -1,
    scale_rows = TRUE,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_cols = c("condition", "cage"),
    show_rownames = TRUE,
    title_text = "Gene Expression Heatmap"
) {
  # 1. Validate expression matrix and metadata
  if (!is.matrix(expression_matrix) ||
      !is.numeric(expression_matrix)) {
    stop("expression_matrix must be a numeric matrix.")
  }

  if (any(!is.finite(expression_matrix))) {
    stop("Expression matrix contains missing or non-finite values.")
  }

  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame.")
  }

  if (!identical(colnames(expression_matrix),
                 rownames(metadata))) {
    stop("Sample names and order must match.")
  }

  if (!is.numeric(top_n) || length(top_n) != 1 ||
      is.na(top_n) || !is.finite(top_n) ||
      top_n < 1 || top_n != floor(top_n)) {
    stop("top_n must be a positive whole number.")
  }


  # 2. Select genes for heatmap

  if (!is.null(gene_list)) {

    selected_genes <- unique(as.character(gene_list))

  } else if (!is.null(de_results)) {

    required_cols <- c(
      gene_id_col,
      gene_name_col,
      "log2FoldChange",
      "pvalue"
    )

    missing_cols <- setdiff(required_cols, names(de_results))

    if (length(missing_cols) > 0) {
      stop(
        "Missing columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }

    de_df <- de_results

    # Convert identifiers to character
    de_df[[gene_id_col]] <- as.character(
      de_df[[gene_id_col]]
    )

    # Keep genes present in the expression matrix
    de_df <- de_df[
      !is.na(de_df[[gene_id_col]]) &
        de_df[[gene_id_col]] %in%
        rownames(expression_matrix),
      ,
      drop = FALSE
    ]

    # Keep valid statistical results
    de_df <- de_df[
      is.finite(de_df$pvalue) &
        is.finite(de_df$log2FoldChange) &
        de_df$pvalue < p_thr,
      ,
      drop = FALSE
    ]

    # Avoid duplicate gene identifiers
    de_df <- de_df[
      !duplicated(de_df[[gene_id_col]]),
      ,
      drop = FALSE
    ]

    # Select upregulated genes
    top_up <- de_df[
      de_df$log2FoldChange >= up_thr,
      ,
      drop = FALSE
    ]

    top_up <- top_up[
      order(
        top_up$pvalue,
        -top_up$log2FoldChange
      ),
      ,
      drop = FALSE
    ]

    top_up <- utils::head(top_up, top_n)

    # Select downregulated genes
    top_down <- de_df[
      de_df$log2FoldChange <= down_thr,
      ,
      drop = FALSE
    ]

    top_down <- top_down[
      order(
        top_down$pvalue,
        top_down$log2FoldChange
      ),
      ,
      drop = FALSE
    ]

    top_down <- utils::head(top_down, top_n)

    # Combine the two gene groups
    selected_df <- rbind(top_up, top_down)

    selected_genes <- as.character(
      selected_df[[gene_id_col]]
    )

  } else {
    stop("Provide either gene_list or de_results.")
  }

  # 3. Match genes to expression matrix

  selected_genes <- selected_genes[
    !is.na(selected_genes) &
      nzchar(selected_genes)
  ]

  selected_genes <- unique(selected_genes)

  selected_genes <- selected_genes[
    selected_genes %in% rownames(expression_matrix)
  ]

  if (length(selected_genes) < 2) {
    stop("At least two matching genes are required.")
  }

  heatmap_matrix <- expression_matrix[
    selected_genes,
    ,
    drop = FALSE
  ]

  # 3a. Convert Entrez IDs to gene symbols

  if (!is.null(de_results)) {

    gene_ids <- as.character(
      de_results[[gene_id_col]]
    )

    gene_symbols <- as.character(
      de_results[[gene_name_col]]
    )

    matched_indices <- match(
      rownames(heatmap_matrix),
      gene_ids
    )

    display_names <- gene_symbols[matched_indices]

    # Use Entrez IDs if symbols are missing
    missing_names <- is.na(display_names) |
      !nzchar(trimws(display_names))

    display_names[missing_names] <-
      rownames(heatmap_matrix)[missing_names]

    # Ensure unique labels
    rownames(heatmap_matrix) <- make.unique(
      display_names
    )
  }


  # 4. Apply row-wise Z-score normalization
  if (scale_rows) {

    gene_sd <- apply(heatmap_matrix, 1, stats::sd)

    heatmap_matrix <- heatmap_matrix[
      is.finite(gene_sd) & gene_sd > 0,
      ,
      drop = FALSE
    ]

    if (nrow(heatmap_matrix) < 2) {
      stop("At least two genes with nonzero variance are required.")
    }

    heatmap_matrix <- t(scale(t(heatmap_matrix)))
  }

  # 5. Prepare sample annotations
  if (is.null(annotation_cols) ||
      length(annotation_cols) == 0) {

    annotation_df <- NULL

  } else {

    missing_annotations <- setdiff(
      annotation_cols,
      colnames(metadata)
    )

    if (length(missing_annotations) > 0) {
      stop(
        "Missing annotation columns: ",
        paste(missing_annotations, collapse = ", ")
      )
    }

    annotation_df <- metadata[
      , annotation_cols, drop = FALSE
    ]

    # Shorten condition labels for display
    if ("condition" %in% colnames(annotation_df)) {
      annotation_df$condition <- as.character(annotation_df$condition)

      annotation_df$condition[
        annotation_df$condition == "SCN5A_KD"
      ] <- "KD"

      annotation_df$condition <- factor(
        annotation_df$condition,
        levels = c("Control", "KD")
      )
    }

    annotation_df[] <- lapply(
      annotation_df,
      function(x) {
        if (is.character(x)) factor(x) else x
      }
    )
  }

  # 6. Define heatmap colors
  heatmap_colors <- grDevices::colorRampPalette(
    c("blue", "white", "red")
  )(100)

  # 7. Generate heatmap

  # Generate shorter sample labels
  sample_labels <- colnames(heatmap_matrix)

  sample_labels <- sub(
    "^ctrl\\.([0-9]+)\\.hg38\\.bam$",
    "Control_\\1",
    sample_labels
  )

  sample_labels <- sub(
    "^sh70c\\.([0-9]+)\\.hg38\\.bam$",
    "KD_\\1",
    sample_labels
  )

  p <- pheatmap::pheatmap(
    mat = heatmap_matrix,
    color = heatmap_colors,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_col = annotation_df,
    show_rownames = show_rownames,
    main = title_text,
    labels_col = sample_labels,
    border_color = NA,
    silent = TRUE,
    legend = TRUE,
    legend_breaks = c(-1.5, -1, 0, 1, 1.5),
    legend_labels = c("-1.5", "-1", "0", "1", "1.5")
  )

  grid::grid.newpage()
  grid::grid.draw(p$gtable)

  return(invisible(p))
  }
