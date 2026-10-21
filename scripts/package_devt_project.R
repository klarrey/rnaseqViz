## Differential gene expression analysis using DESeq2



# Install GEOquery if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}

library(GEOquery)

# Select dataset
gse_id <- "GSE228621"

# Retrieve GEO Series metadata
gse <- getGEO(
  gse_id,
  GSEMatrix = FALSE
)

# View study information
Meta(gse)$title
Meta(gse)$summary
Meta(gse)$overall_design

# View sample identifiers
names(GSMList(gse))

# View supplementary file information
Meta(gse)$supplementary_file


library(GEOquery)

# Create a directory for our project data
dir.create(
  "data",
  showWarnings = FALSE
)

# Download supplementary files
getGEOSuppFiles(
  GEO = "GSE228621",
  baseDir = "data"
)


# Define file path
count_file <- file.path(
  "data",
  "GSE228621",
  "GSE228621_CountMatrix.csv.gz"
)

# Read compressed CSV
counts <- read.csv(
  gzfile(count_file),
  row.names = 1,
  check.names = FALSE
)

# Inspect dimensions
dim(counts)

# Inspect first few rows
head(counts)

# Inspect sample names
colnames(counts)

# Inspect data types
str(counts)



# Extract individual GEO sample records
gsm_list <- GSMList(gse)

# Create a sample metadata table
sample_info <- data.frame(
  GSM = names(gsm_list),

  title = sapply(
    gsm_list,
    function(x) Meta(x)$title
  ),

  characteristics = sapply(
    gsm_list,
    function(x) paste(
      Meta(x)$characteristics_ch1,
      collapse = "; "
    )
  ),

  row.names = NULL
)

# Display metadata
print(sample_info)



# Extract gene annotations
gene_annotations <- counts[, c("gene_name", "entrez_id")]

# Extract only RNA-seq count columns
count_matrix <- counts[, 1:12]

# Inspect dimensions
dim(count_matrix)

# Inspect first few rows
head(count_matrix)



# Convert data frame to matrix
count_matrix <- as.matrix(count_matrix)

# Check data type
typeof(count_matrix)

# Check for missing values
sum(is.na(count_matrix))

# Check for negative values
sum(count_matrix < 0)

# Check whether all values are integers
all(count_matrix == round(count_matrix))


# Create sample metadata
coldata <- data.frame(
  sample = colnames(count_matrix),

  condition = c(
    rep("Control", 6),
    rep("SCN5A_KD", 6)
  ),

  row.names = colnames(count_matrix)
)

# Convert condition to factor
coldata$condition <- factor(
  coldata$condition,
  levels = c("Control", "SCN5A_KD")
)

# Display metadata
print(coldata)


# Verify sample alignment
all(
  colnames(count_matrix) == rownames(coldata)
)

## This check is essential because DESeq2 requires the count matrix
## columns and metadata rows to correspond to the same samples in the same order.



# Extract cage identifiers from GEO metadata
coldata$cage <- sub(
  ".*cage: ([^;]+).*",
  "\\1",
  sample_info$characteristics
)

# Display sample metadata
print(coldata)

# Count samples per cage and condition
table(coldata$cage, coldata$condition)




# 1. Check the count matrix
dim(count_matrix)
str(count_matrix)

# 2. Verify count integrity
sum(is.na(count_matrix))
sum(count_matrix < 0)
all(count_matrix == round(count_matrix))

# 3. Examine cage and treatment distribution
table(coldata$cage, coldata$condition)


##These results indicate that the count matrix has
##the expected structure and contains valid nonnegative integer values. Note: DESeq2 does not
##accept negative and missing values

# Check sample order
stopifnot(
  identical(
    colnames(count_matrix),
    rownames(coldata)
  )
)

# Check treatment distribution
table(coldata$condition)

# Check for duplicated gene identifiers
sum(duplicated(rownames(count_matrix)))

# Check for duplicated sample identifiers
anyDuplicated(colnames(count_matrix))



#Install DESeq2
BiocManager::install("DESeq2")

#creat DESeq2 object
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ condition
)

dds



# Identify genes meeting the filtering criterion
keep <- rowSums(counts(dds) >= 10) >= 3

# Number of genes retained
sum(keep)

# Number of genes removed
sum(!keep)

# Apply filtering
dds <- dds[keep, ]

# Inspect filtered dataset
dim(dds)



# Perform differential expression analysis
dds <- DESeq(dds)

# Extract SCN5A knockdown vs control results
res <- results(
  dds,
  contrast = c(
    "condition",
    "SCN5A_KD",
    "Control"
  )
)

# Summarize results
summary(res)

# Convert results to a data frame
res_df <- as.data.frame(res)

# Add gene identifiers
res_df$entrez_id <- rownames(res_df)

# Inspect results
head(res_df)


#create a folder call results
dir.create("results", showWarnings = FALSE)

# Save differential expression results
write.csv(
  res_df,
  "results/GSE228621_DESeq2_results.csv",
  row.names = FALSE
)

# Save the DESeq2 object
saveRDS(
  dds,
  "results/GSE228621_dds.rds"
)


summary(res)
dim(dds)



# Define thresholds outside the function
up_thr   <- 1
down_thr <- -1
p_thr    <- 0.05

# Convert DESeq2 results to data frame
#res_df <- as.data.frame(res)

# Count upregulated genes
up_n <- sum(
  res_df$log2FoldChange >= up_thr &
    res_df$pvalue < p_thr,
  na.rm = TRUE
)

# Count downregulated genes
down_n <- sum(
  res_df$log2FoldChange <= down_thr &
    res_df$pvalue < p_thr,
  na.rm = TRUE
)

cat("Upregulated:", up_n, "\n")
cat("Downregulated:", down_n, "\n")


# Next step: quality-control visualizations
#library(DESeq2)
library(ggplot2)

# Apply variance-stabilizing transformation
vsd <- vst(dds, blind = TRUE)

# Generate PCA plot
pca_plot <- plotPCA(
  vsd,
  intgroup = "condition"
)

print(pca_plot)




#Next step: PCA colored by cage

# Extract PCA coordinates
pca_data <- plotPCA(
  vsd,
  intgroup = c("condition", "cage"),
  returnData = TRUE
)

# Extract variance explained
percentVar <- round(
  100 * attr(pca_data, "percentVar")
)

# Generate PCA plot
pca_cage <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = cage,
    shape = condition
  )
) +
  geom_point(size = 4) +
  labs(
    title = "PCA of Breast Cancer RNA-seq Samples",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    color = "Cage",
    shape = "Condition"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

print(pca_cage)



library(DESeq2)

# Model 1: Condition only
dds_condition <- dds

design(dds_condition) <- ~ condition

dds_condition <- DESeq(dds_condition)

res_condition <- results(
  dds_condition,
  contrast = c("condition", "SCN5A_KD", "Control")
)

# Model 2: Cage-adjusted
dds_cage <- dds

colData(dds_cage)$cage <- factor(
  colData(dds_cage)$cage
)

design(dds_cage) <- ~ cage + condition

dds_cage <- DESeq(dds_cage)

res_cage <- results(
  dds_cage,
  contrast = c("condition", "SCN5A_KD", "Control")
)




res_condition <- res
# Compare results
summary(res_condition)
summary(res_cage)



# Compare significant gene counts for the two models

# FDR threshold
padj_thr <- 0.05

# Significant genes in condition-only model
sig_condition <- sum(
  res_condition$padj < padj_thr,
  na.rm = TRUE
)

# Significant genes in cage-adjusted model
sig_cage <- sum(
  res_cage$padj < padj_thr,
  na.rm = TRUE
)

data.frame(
  Model = c("Condition only", "Cage adjusted"),
  Significant_genes = c(sig_condition, sig_cage)
)



## Compare the estimated log2 fold changes

library(dplyr)
library(ggplot2)
library(tidyverse)

# Create comparison table
comparison_df <- data.frame(
  entrez_id = rownames(res_condition),
  logFC_condition = res_condition$log2FoldChange,
  logFC_cage = res_cage$log2FoldChange,
  padj_condition = res_condition$padj,
  padj_cage = res_cage$padj
)

# Remove genes with missing fold changes
comparison_df <- comparison_df %>%
  filter(
    is.finite(logFC_condition),
    is.finite(logFC_cage)
  )

# Plot comparison
ggplot(
  comparison_df,
  aes(
    x = logFC_condition,
    y = logFC_cage
  )
) +
  geom_point(
    alpha = 0.4,
    size = 1
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "red"
  ) +
  labs(
    title = "Comparison of DESeq2 Model Estimates",
    x = "log2FC: Condition-only model",
    y = "log2FC: Cage-adjusted model"
  ) +
  theme_bw(base_size = 14)


cor(
  comparison_df$logFC_condition,
  comparison_df$logFC_cage,
  method = "pearson"
)




library(dplyr)
library(ggplot2)

# Create comparison table
se_comparison <- data.frame(
  entrez_id = rownames(res_condition),
  SE_condition = res_condition$lfcSE,
  SE_cage = res_cage$lfcSE
)

# Remove missing values
se_comparison <- se_comparison %>%
  filter(
    is.finite(SE_condition),
    is.finite(SE_cage)
  )

# Compare standard errors
ggplot(
  se_comparison,
  aes(x = SE_condition, y = SE_cage)
) +
  geom_point(alpha = 0.3, size = 1) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "red"
  ) +
  labs(
    title = "Comparison of DESeq2 Standard Errors",
    x = "Standard error: Condition-only model",
    y = "Standard error: Cage-adjusted model"
  ) +
  theme_bw(base_size = 14)


# Quantify the differences
# Median standard errors
median(se_comparison$SE_condition)
median(se_comparison$SE_cage)

# Percentage of genes with increased standard errors
mean(
  se_comparison$SE_cage >
    se_comparison$SE_condition
) * 100



#Check whether the models identify the same genes
# Genes significant in condition-only model
sig_genes_condition <- rownames(res_condition)[
  !is.na(res_condition$padj) &
    res_condition$padj < 0.05
]

# Genes significant in cage-adjusted model
sig_genes_cage <- rownames(res_cage)[
  !is.na(res_cage$padj) &
    res_cage$padj < 0.05
]

# Shared significant genes
shared_genes <- intersect(
  sig_genes_condition,
  sig_genes_cage
)

# Count genes
length(shared_genes)

# Unique to each model
length(setdiff(sig_genes_condition, sig_genes_cage))
length(setdiff(sig_genes_cage, sig_genes_condition))




# Identify genes unique to the cage-adjusted model
unique_cage_genes <- setdiff(
  sig_genes_cage,
  sig_genes_condition
)

# Extract their statistics from both models
unique_cage_results <- data.frame(
  entrez_id = unique_cage_genes,

  logFC_condition =
    res_condition[unique_cage_genes, "log2FoldChange"],

  padj_condition =
    res_condition[unique_cage_genes, "padj"],

  logFC_cage =
    res_cage[unique_cage_genes, "log2FoldChange"],

  padj_cage =
    res_cage[unique_cage_genes, "padj"]
)

print(unique_cage_results)





###Here I added the gene symbol from the original count matrix data to the analysis results dataframe

# Extract gene annotations from original data
gene_annotations <- data.frame(
  entrez_id = as.character(counts$entrez_id),
  gene_name = as.character(counts$gene_name)
)

# Ensure DESeq2 results contain Entrez IDs
res_df$entrez_id <- rownames(res_df)

# Match gene symbols using Entrez IDs
res_df$gene_name <- gene_annotations$gene_name[
  match(
    res_df$entrez_id,
    gene_annotations$entrez_id
  )
]

# Inspect the results
head(res_df[, c(
  "entrez_id",
  "gene_name",
  "log2FoldChange",
  "pvalue",
  "padj"
)])



# Count genes without symbols
sum(is.na(res_df$gene_name))

# Check first 10 genes
head(
  res_df[, c("entrez_id", "gene_name")],
  10
)

# Check for duplicate Entrez IDs
sum(duplicated(gene_annotations$entrez_id))


saveRDS(
  res_df,
  "results/GSE228621_DESeq2_annotated_results.rds"
)

### function to create volcano plot

library(ggrepel)

make_volcano <- function(df, title_text,
                         up_thr, down_thr, p_thr,
                         top_n = 10) {

  # Validate thresholds
  stopifnot(
    up_thr > 0,
    down_thr < 0,
    p_thr > 0,
    p_thr < 1
  )

  # Prepare differential expression results
  df <- res_df %>%
    mutate(
      log2FoldChange = as.numeric(log2FoldChange),
      pvalue = as.numeric(pvalue),

      negLogP = -log10(pvalue),

      status = case_when(
        log2FoldChange >= up_thr & pvalue < p_thr ~ "Up",
        log2FoldChange <= down_thr & pvalue < p_thr ~ "Down",
        TRUE ~ "NS"
      )
    )

  # Count genes in each category
  total_n <- nrow(df)

  up_n <- sum(df$status == "Up", na.rm = TRUE)

  down_n <- sum(df$status == "Down", na.rm = TRUE)


  #Exclude genes without valid symbols from labeling
  # Select top 10 upregulated genes
  top_up <- df %>%
    filter(
      status == "Up",
      !is.na(gene_name),
      gene_name != ""
    ) %>%
    arrange(pvalue, desc(log2FoldChange)) %>%
    slice_head(n = top_n)


  # Select top 10 downregulated genes
  top_down <- df %>%
    filter(
      status == "Down",
      !is.na(gene_name),
      gene_name != ""
    ) %>%
    arrange(pvalue, log2FoldChange) %>%
    slice_head(n = top_n)

  # Combine genes for labeling
  label_df <- bind_rows(top_up, top_down)

  # Determine symmetric x-axis limits
  max_abs <- max(
    abs(df$log2FoldChange),
    na.rm = TRUE
  )

  # Generate volcano plot
  g <- ggplot(
    df,
    aes(x = log2FoldChange, y = negLogP)
  ) +

    geom_point(
      aes(color = status),
      alpha = 0.7,
      size = 1.8
    ) +

    scale_color_manual(
      values = c(
        "Up" = "red",
        "Down" = "blue",
        "NS" = "grey70"
      ),
      breaks = c("Up", "Down", "NS")
    ) +

    geom_vline(
      xintercept = c(down_thr, up_thr),
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.5
    ) +

    geom_hline(
      yintercept = -log10(p_thr),
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.5
    ) +

    geom_text_repel(
      data = label_df,
      aes(label = gene_name),
      size = 3.5,
      max.overlaps = Inf,
      box.padding = 0.8,
      point.padding = 0.5,
      force = 3,
      force_pull = 0.5,
      max.iter = 20000,
      min.segment.length = 0,
      segment.color = "grey50",
      seed = 123
    ) +

    labs(
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

    coord_cartesian(
      xlim = c(-max_abs, max_abs)
    ) +

    theme_bw(base_size = 15) +

    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )

  return(g)
}

# Thresholds
up_thr   <- 1
down_thr <- -1
p_thr    <- 0.055

# Number of genes to label in each category
top_n <- 10



# Read differential expression results
res_df <- read.csv("GSE228621_DESeq2_results.csv")

# Generate plot
volcano_plot <- make_volcano(
  df = res_df,
  title_text = "SCN5A Knockdown vs Control",
  up_thr = up_thr,
  down_thr = down_thr,
  p_thr = p_thr,
  top_n = top_n
)

# Display plot
print(volcano_plot)

#save the plot

ggsave(
  filename = "SCN5A_KD_volcano.png",
  plot = volcano_plot,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)





### Develop the rnaseqViz package

install.packages(c(
  "usethis",
  "devtools",
  "roxygen2",
  "testthat"
))



library(usethis)

create_package(
  path = "C:/Users/larre/OneDrive/Desktop/portfolio_projects/rnaseqViz",
  open = TRUE
)



