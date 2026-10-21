test_that("plot_volcano returns a ggplot object", {
  # Small example differential expression dataset
  test_df <- data.frame(
    gene_name = c("GENE1", "GENE2", "GENE3",
                  "GENE4", "GENE5"),

    log2FoldChange = c(2, -2, 0.5, 1.5, -1.5),

    pvalue = c(0.001, 0.002, 0.5, 0.01, 0.03)
  )

  # Generate volcano plot
  p <- plot_volcano(
    df = test_df,
    title_text = "Test Volcano Plot",
    up_thr = 1,
    down_thr = -1,
    p_thr = 0.05,
    top_n = 2
  )

  # Verify output
  expect_s3_class(p, "ggplot")
})


#test for gene classification
test_that("plot_volcano classifies genes correctly", {

  test_df <- data.frame(
    gene_name = c("GENE1", "GENE2", "GENE3",
                  "GENE4", "GENE5"),
    log2FoldChange = c(2, -2, 0.5, 1.5, -1.5),
    pvalue = c(0.001, 0.002, 0.5, 0.01, 0.03)
  )

  p <- plot_volcano(
    df = test_df,
    title_text = "Test Volcano Plot",
    up_thr = 1,
    down_thr = -1,
    p_thr = 0.05,
    top_n = 2
  )

  # Check classifications stored in the plot data
  expect_equal(
    sum(p$data$status == "Up"),
    2
  )

  expect_equal(
    sum(p$data$status == "Down"),
    2
  )

  expect_equal(
    sum(p$data$status == "NS"),
    1
  )
})


#test for adjustable thresholds
test_that("plot_volcano supports asymmetric thresholds", {

  test_df <- data.frame(
    gene_name = c("A", "B", "C", "D"),
    log2FoldChange = c(1.5, 2.5, -1.5, -2.5),
    pvalue = c(0.01, 0.01, 0.01, 0.01)
  )

  p <- plot_volcano(
    df = test_df,
    title_text = "Asymmetric Threshold Test",
    up_thr = 2,
    down_thr = -1,
    p_thr = 0.05,
    top_n = 2
  )

  expect_equal(sum(p$data$status == "Up"), 1)
  expect_equal(sum(p$data$status == "Down"), 2)
  expect_equal(sum(p$data$status == "NS"), 1)
})



test_that("plot_volcano rejects missing columns", {

  bad_df <- data.frame(
    gene_name = c("A", "B"),
    log2FoldChange = c(2, -2)
  )

  expect_error(
    plot_volcano(bad_df, "Test"),
    "Missing required columns: pvalue"
  )
})

test_that("plot_volcano rejects invalid thresholds", {

  test_df <- data.frame(
    gene_name = c("A", "B"),
    log2FoldChange = c(2, -2),
    pvalue = c(0.01, 0.02)
  )

  expect_error(
    plot_volcano(test_df, "Test", up_thr = -1),
    "up_thr must be a positive number"
  )

  expect_error(
    plot_volcano(test_df, "Test", down_thr = 1),
    "down_thr must be a negative number"
  )

  expect_error(
    plot_volcano(test_df, "Test", p_thr = 1.5),
    "p_thr must be between 0 and 1"
  )
})
