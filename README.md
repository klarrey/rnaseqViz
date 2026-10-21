# rnaseqViz

### An R Package and Interactive Shiny Application for RNA-seq Data Visualization

`rnaseqViz` is an R package that provides reusable functions and an interactive Shiny application for exploring and visualizing bulk RNA-seq differential gene expression results.

The application currently uses a breast cancer RNA-seq dataset from NCBI Gene Expression Omnibus (GEO), accession **GSE228621**, investigating gene expression changes following SCN5A knockdown in MDA-MB-231 breast cancer xenograft tumors.

The project demonstrates an end-to-end bioinformatics visualization workflow integrating DESeq2 results, exploratory data analysis, interactive visualizations, and downloadable results.

## Features

**1. Principal Component Analysis (PCA)**

- Explore variation and clustering among RNA-seq samples.
- Color samples by experimental condition.
- Differentiate samples by cage using point shapes.
- Interactively select the number of variable genes.
- Export PCA visualizations as PNG and PDF.

**2. Volcano Plot**

- Visualize upregulated, downregulated, and nonsignificant genes.
- Customize log2 fold-change and p-value thresholds.
- Adjust the number of labeled genes.
- Display dynamic differential expression counts.
- Export publication-quality plots.

**3. Gene Expression Heatmap**

- Visualize the most upregulated and downregulated genes.
- Adjust differential expression thresholds and the number of selected genes.
- Display row-scaled gene expression patterns.
- Annotate samples by experimental condition and cage.
- Export heatmaps as PNG and PDF.

**4. Differential Expression Results**

- Explore DESeq2 results in a searchable, sortable table.
- Filter genes by expression thresholds and classification.
- View dynamic summaries of upregulated, downregulated, and nonsignificant genes.
- Download filtered differential expression results as CSV.

**5. Study Documentation**

- Review the dataset, experimental design, and analysis methodology.
- Access descriptions of the application's visualization and analysis features.

## Dataset

| Property | Description |
|---|---|
| GEO accession | GSE228621 |
| Biological model | MDA-MB-231 breast cancer xenograft tumors |
| Experimental comparison | SCN5A knockdown vs Control |
| Total samples | 12 |
| Control samples | 6 |
| SCN5A knockdown samples | 6 |
| Genes analyzed | 14,703 |
| Differential expression analysis | DESeq2 |

At the default thresholds of absolute log2 fold change ≥ 1 and unadjusted p-value < 0.05, the analysis identifies 32 upregulated and 110 downregulated genes.

These counts reflect the stated thresholds and should not be interpreted as FDR-adjusted discoveries.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")

remotes::install_github(
  "klarrey/rnaseqViz",
  dependencies = TRUE
)
```

## Launch the Shiny Application

```r
library(rnaseqViz)

rnaseqViz::run_app()
```

The application provides interactive access to the included breast cancer RNA-seq dataset.

## Using the R Functions

The package also provides reusable visualization functions.

For example:

```r
library(rnaseqViz)

# View available exported functions
getNamespaceExports("rnaseqViz")

# Review PCA function arguments
args(rnaseqViz::plot_pca)

# Review volcano plot function arguments
args(rnaseqViz::plot_volcano)

# Review heatmap function arguments
args(rnaseqViz::plot_heatmap)
```

The plotting functions can be integrated into other R-based bioinformatics analysis workflows.

## Technologies

- R
- Shiny
- bslib
- ggplot2
- DESeq2
- DT
- pheatmap
- testthat
- roxygen2

## Project Status

**Development version: 0.0.0.9000**

The first application version focuses on a single breast cancer RNA-seq dataset. Future development may include user-uploaded datasets, additional analysis modules, and broader support for transcriptomics workflows.

## Author

Enoch Larrey, Ph.D.

Bioinformatics | Computational Biology | RNA-seq Data Analysis | R Package Development

## Data Source

NCBI Gene Expression Omnibus: [GSE228621](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228621)
