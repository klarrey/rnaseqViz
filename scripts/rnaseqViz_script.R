
#set the working directory
setwd("C:/Users/larre/OneDrive/Desktop/portfolio_projects/rnaseqViz")


## Start creating the package files

# Check current working directory
getwd()

# List package files
list.files()

# Confirm DESCRIPTION exists
file.exists("DESCRIPTION")

# Confirm R folder exists
dir.exists("R")



# Generate the MIT license
usethis::use_mit_license(
  copyright_holder = "Enoch Larrey"
)

##Add the required dependencies
usethis::use_package("ggplot2", type = "Imports")
usethis::use_package("dplyr", type = "Imports")
usethis::use_package("ggrepel", type = "Imports")


# Display DESCRIPTION contents
read.dcf("DESCRIPTION")

# Check that package metadata can be read
packageDescription("rnaseqViz", lib.loc = NULL)


#Create a new R script
usethis::use_r("plot_volcano")

#Generate documentation
devtools::document()


devtools::load_all()

#Check whether the function is available
exists("plot_volcano")

#open the function's documentation
?plot_volcano



#Test plot_volcano() using our breast cancer dataset

## Load the DESeq2 results
res_df <- readRDS("C:/Users/larre/OneDrive/Desktop/portfolio_projects/results/GSE228621_DESeq2_annotated_results.rds")

colnames(res_df)


#Verify the required columns
head(res_df[, c(
  "entrez_id",
  "gene_name",
  "log2FoldChange",
  "pvalue",
  "padj"
)])

#Test the package function
# Load package in development mode
devtools::load_all()

# Generate volcano plot
volcano_plot <- plot_volcano(
  df = res_df,
  title_text = "SCN5A Knockdown vs Control",
  up_thr = 1,
  down_thr = -1,
  p_thr = 0.055,
  top_n = 10
)

# Display plot
print(volcano_plot)

#Check whether the plot object was created
class(volcano_plot)

inherits(volcano_plot, "ggplot")

#Initialize the testing framework
usethis::use_testthat(edition = 3)


#Create a test file
usethis::use_test("plot_volcano")


#Run the tests
devtools::test()


##Add input validation


#update the package documentation and NAMESPACE based on current roxygen2 comments
devtools::document()


#Run the package check
devtools::check()

# Test the modifications
devtools::document()
devtools::test()
devtools::check()


#Resolve the scripts directory note
usethis::use_build_ignore("scripts")

readLines(".Rbuildignore")


#Resolve the undefined global variable note
usethis::use_package("rlang", type = "Imports")




#============================================================
#Building pca package
#============================================================

#load library
library(DESeq2)

exists("dds")


## Load the DESeqDataSet object
dds <- readRDS("C:/Users/larre/OneDrive/Desktop/portfolio_projects/results/GSE228621_dds.rds")
head(dds$cage)


#Perform variance-stabilizing transformation
vsd <- vst(dds, blind = FALSE)

#extract the expression matrix
expression_matrix <- assay(vsd)


#Prepare sample metadata
sample_metadata <- as.data.frame(
  colData(vsd)
)


#verify sample alignment
identical(
  colnames(expression_matrix),
  rownames(sample_metadata)
)


#Create the R package function file
usethis::use_r("plot_pca")


head(sample_metadata)
colnames(sample_metadata)

#confirm that ggplot2 and rlang are declared as dependencies
usethis::use_package("ggplot2", type = "Imports")
usethis::use_package("rlang", type = "Imports")


#Generate documentation and test with breast cancer dataset
devtools::document()
devtools::load_all()


#generate the PCA plot
pca_plot <- plot_pca(
  expression_matrix = expression_matrix,
  metadata = sample_metadata,
  color_by = "condition",
  shape_by = "cage",
  top_n = 500,
  title_text = "PCA of Breast Cancer RNA-seq Samples"
)

print(pca_plot)


#Create automated tests for plot_pca()
usethis::use_test("plot_pca")



devtools::document()
devtools::test()
devtools::check()





#######################################################
#Design the heatmap function
#######################################################
usethis::use_package("pheatmap", type = "Imports")


usethis::use_r("plot_heatmap")


#load library
library(DESeq2)



## Load the DESeqDataSet object
dds <- readRDS("C:/Users/larre/OneDrive/Desktop/portfolio_projects/results/GSE228621_dds.rds")

#Perform variance-stabilizing transformation
vsd <- vst(dds, blind = FALSE)

#extract the expression matrix
expression_matrix <- assay(vsd)


#Prepare sample metadata
sample_metadata <- as.data.frame(
  colData(vsd)
)


#verify sample alignment
identical(
  colnames(expression_matrix),
  rownames(sample_metadata)
)



## Load the DESeq2 results
res_df <- readRDS("C:/Users/larre/OneDrive/Desktop/portfolio_projects/results/GSE228621_DESeq2_annotated_results.rds")

#

# Define file path
count_file <- "C:/Users/larre/OneDrive/Desktop/portfolio_projects/data/GSE228621/GSE228621_CountMatrix.csv.gz"

# Read compressed CSV
counts <- read.csv(
  gzfile(count_file),
  row.names = 1,
  check.names = FALSE
)

#verify the presence and matching of gene_name and entrez_id columns
colnames(res_df)
head(rownames(expression_matrix))
head(counts[, c("gene_name", "entrez_id")])


# Number of genes in expression matrix
nrow(expression_matrix)

# Number of genes in DESeq2 results
nrow(res_df)

# Count matching Entrez IDs
sum(
  as.character(res_df$entrez_id) %in%
    rownames(expression_matrix)
)

# Check for duplicated Entrez IDs
sum(duplicated(as.character(res_df$entrez_id)))

# Check for missing gene symbols
sum(
  is.na(res_df$gene_name) |
    trimws(as.character(res_df$gene_name)) == ""
)



devtools::document()
devtools::load_all()

heatmap_plot <- plot_heatmap(
  expression_matrix = expression_matrix,
  metadata = sample_metadata,
  de_results = res_df,
  gene_id_col = "entrez_id",
  gene_name_col = "gene_name",
  top_n = 10,
  p_thr = 0.055,
  up_thr = 1,
  down_thr = -1,
  scale_rows = TRUE,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_cols = c("condition", "cage"),
  show_rownames = TRUE,
  title_text = "SCN5A Knockdown vs Control"
)



code <- readLines("R/plot_heatmap.R")

cat(
  paste(
    seq_len(min(55, length(code))),
    utils::head(code, 55),
    sep = ": ",
    collapse = "\n"
  )
)



# Number of upregulated genes
sum(
  res_df$log2FoldChange >= 1 &
    res_df$pvalue < 0.055,
  na.rm = TRUE
)

# Number of downregulated genes
sum(
  res_df$log2FoldChange <= -1 &
    res_df$pvalue < 0.055,
  na.rm = TRUE
)



res_df[
  as.character(res_df$entrez_id) %in%
    c("286499", "79047", "4110", "4066", "972"),
  c("entrez_id", "gene_name", "log2FoldChange", "pvalue")
]


# Check that 20 genes were selected
nrow(heatmap_plot$tree_row$labels)

# Display the gene symbols
heatmap_plot$tree_row$labels


#automated testing of package
usethis::use_test("plot_heatmap")


devtools::test()
devtools::document()
devtools::check()



#######################################################################
#Building the shiny app
#######################################################################


#Create the Shiny application structure
usethis::use_directory("inst/shiny")

# Verify the required objects exist
exists("res_df")
exists("expression_matrix")
exists("sample_metadata")


colnames(res_df)



# Save datasets as internal package data
usethis::use_data(
  res_df,
  expression_matrix,
  sample_metadata,
  internal = TRUE,
  overwrite = TRUE
)



#Create the Shiny application structure
#Install the required packages
install.packages(c(
  "shiny",
  "bslib",
  "DT"
))


#Add dependencies to your R package
usethis::use_package("shiny")
usethis::use_package("bslib")
usethis::use_package("DT")


#Create the Shiny application directory
dir.create(
  "inst/shiny",
  recursive = TRUE,
  showWarnings = FALSE
)


#verify that your packaged data can be loaded
devtools::load_all()

dim(res_df)
dim(expression_matrix)
dim(sample_metadata)




#Create the first working Shiny dashboard
#Create app.R
file.create("inst/shiny/app.R")


#Load the package and run the app
devtools::load_all()

shiny::runApp("inst/shiny")


#Add the PCA tab
args(plot_pca)


class(expression_matrix)


#============================================================
#Add Differential Expression Results tab
#============================================================

#Step 1. Install the DT package
install.packages("DT")

#open the rnaseqViz package project and register DT as a dependency
usethis::use_package("DT")

#Step 2. Add the DE Results navigation tab

#==========================================================
#interactive differential expression filtering
#==========================================================
#Step 1. Add filtering controls to the DE Results tab

#Test the application
devtools::load_all()

shiny::runApp("inst/shiny")



#=========================================================
#Step 1 — Validate the R package
#=========================================================

#1.1 Open your package project

#1.2 Run the complete validation
# Load the current package
devtools::load_all()

# Update package documentation
devtools::document()

# Run automated tests
devtools::test()

# Run R CMD check
devtools::check()


#==================================================================
#Create a Shiny application launcher
#==================================================================

usethis::use_r("run_app")

#pasted some code to run_app.R file


#launch the application after installation
library(rnaseqViz)

run_app()


#==============================================================
#Step:2 Prepare the GitHub repository
#==============================================================
#2.1 Create README.md

usethis::use_readme_md()


#2.2 Initialize Git
usethis::use_git()


#Step 3.2 — Check the Git repository
getwd()

system("git status")


#Step 3.3 — Review your .gitignore
usethis::use_git_ignore(
  c(
    ".Rhistory",
    ".RData",
    ".Ruserdata",
    ".Rproj.user",
    ".Rcheck",
    ".Renviron",
    ".Rprofile"
  )
)


#Step 3.4 — Verify the files before committing
# Files Git is not ignoring
system("git status --short")

# Files Git is configured to ignore
system("git status --ignored --short")



# Remove generated figures and CSVs from Git's index
# Keep the actual files on your computer
system("git rm -r --cached -- scripts/figs")

# Stage the updated .gitignore
system("git add .gitignore")

# Verify which files are staged
system("git status --short")

#Step 3.4 — Create first commit









