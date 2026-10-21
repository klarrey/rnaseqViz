# Create reproducible test data
set.seed(123)

test_matrix <- matrix(
  rnorm(100),
  nrow = 20,
  ncol = 5
)

rownames(test_matrix) <- paste0("Gene", 1:20)
colnames(test_matrix) <- paste0("Sample", 1:5)

test_metadata <- data.frame(
  condition = c("Control", "Control",
                "Treatment", "Treatment", "Treatment"),
  cage = c("A", "B", "A", "B", "C"),
  row.names = colnames(test_matrix)
)

test_that("plot_pca returns a ggplot object", {

  p <- plot_pca(
    expression_matrix = test_matrix,
    metadata = test_metadata,
    color_by = "condition",
    shape_by = "cage",
    top_n = 10
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_pca calculates PC1 and PC2", {

  p <- plot_pca(
    expression_matrix = test_matrix,
    metadata = test_metadata,
    top_n = 10
  )

  expect_true("PC1" %in% names(p$data))
  expect_true("PC2" %in% names(p$data))
  expect_equal(nrow(p$data), 5)
})

test_that("plot_pca supports color and shape variables", {

  p <- plot_pca(
    expression_matrix = test_matrix,
    metadata = test_metadata,
    color_by = "condition",
    shape_by = "cage"
  )

  expect_equal(p$labels$colour, "condition")
  expect_equal(p$labels$shape, "cage")
})

test_that("plot_pca rejects mismatched sample identifiers", {

  bad_metadata <- test_metadata
  rownames(bad_metadata) <- paste0("Wrong", 1:5)

  expect_error(
    plot_pca(test_matrix, bad_metadata),
    "Sample names and order must match"
  )
})

test_that("plot_pca rejects invalid metadata columns", {

  expect_error(
    plot_pca(
      test_matrix,
      test_metadata,
      color_by = "unknown"
    ),
    "color_by column not found"
  )

  expect_error(
    plot_pca(
      test_matrix,
      test_metadata,
      shape_by = "unknown"
    ),
    "shape_by column not found"
  )
})

test_that("plot_pca rejects invalid gene selection", {

  expect_error(
    plot_pca(
      test_matrix,
      test_metadata,
      top_n = -5
    ),
    "top_n must be a whole number"
  )
})

