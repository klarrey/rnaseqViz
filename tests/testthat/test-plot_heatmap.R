
# Create a small synthetic expression matrix
test_expression <- matrix(
  c(
    8, 9, 10, 15, 16, 17,
    20, 19, 18, 10, 9, 8,
    5, 6, 7, 12, 13, 14,
    15, 14, 13, 6, 5, 4
  ),
  nrow = 4,
  byrow = TRUE
)

rownames(test_expression) <- c("101", "102", "103", "104")

colnames(test_expression) <- paste0("Sample", 1:6)

# Create sample metadata
test_metadata <- data.frame(
  condition = rep(c("Control", "KD"), each = 3),
  cage = rep(c("A", "B", "C"), 2),
  row.names = colnames(test_expression)
)

# Create synthetic DESeq2 results
test_de <- data.frame(
  entrez_id = c("101", "102", "103", "104"),
  gene_name = c("GENEA", "GENEB", "GENEC", "GENED"),
  log2FoldChange = c(2, -2, 1.5, -1.5),
  pvalue = c(0.001, 0.002, 0.01, 0.02)
)

test_that("plot_heatmap returns a pheatmap object", {

  p <- plot_heatmap(
    expression_matrix = test_expression,
    metadata = test_metadata,
    de_results = test_de,
    top_n = 2,
    p_thr = 0.05
  )

  expect_s3_class(p, "pheatmap")
})

test_that("plot_heatmap selects genes in both directions", {

  p <- plot_heatmap(
    expression_matrix = test_expression,
    metadata = test_metadata,
    de_results = test_de,
    top_n = 2,
    p_thr = 0.05
  )

  expect_length(p$tree_row$labels, 4)
  expect_setequal(
    p$tree_row$labels,
    c("GENEA", "GENEB", "GENEC", "GENED")
  )
})

test_that("plot_heatmap validates sample alignment", {

  incorrect_metadata <- test_metadata
  rownames(incorrect_metadata) <- paste0("Wrong", 1:6)

  expect_error(
    plot_heatmap(
      expression_matrix = test_expression,
      metadata = incorrect_metadata,
      de_results = test_de
    ),
    "Sample names and order must match"
  )
})
