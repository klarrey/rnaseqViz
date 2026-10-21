
# Load required packages
library(shiny)
library(bslib)



# Note: Internal package data should instead be accessed
# through the namespace. See below.
# Access internal package datasets
res_df <- get("res_df", envir = asNamespace("rnaseqViz"))

expression_matrix <- get(
  "expression_matrix",
  envir = asNamespace("rnaseqViz")
)

sample_metadata <- get(
  "sample_metadata",
  envir = asNamespace("rnaseqViz")
)

# Define user interface
ui <- bslib::page_navbar(
  title = "rnaseqViz",
  theme = bs_theme(bootswatch = "flatly"),

  # Overview tab
  bslib::nav_panel(
    title = "Overview",

    h2("Breast Cancer RNA-seq Analysis"),
    p("SCN5A Knockdown vs Control"),
    p("Dataset: GEO GSE228621"),

    # Summary statistics
    bslib::layout_columns(
      bslib::value_box(
        title = "Genes Analyzed",
        value = format(nrow(res_df), big.mark = ",")
      ),

      bslib::value_box(
        title = "Samples",
        value = "12"
      ),

      bslib::value_box(
        title = "Experimental Groups",
        value = "2"
      )
    ),


    hr(),

    h3("Study Overview"),

    p(
      "This application visualizes differential gene expression",
      "in MDA-MB-231 breast cancer xenograft tumors",
      "following SCN5A knockdown."
    ),

    h3("Experimental Design"),

    tags$ul(
      tags$li("Control: 6 samples"),
      tags$li("SCN5A Knockdown: 6 samples"),
      tags$li("Analysis: DESeq2")
    )
  ),

  bslib::nav_panel(
    title = "PCA Analysis",

    h2("Principal Component Analysis"),

    p(
      "Explore sample variation and clustering",
      "in breast cancer RNA-seq samples."
    ),


    shiny::sliderInput(
      inputId = "top_n",
      label = "Number of variable genes",
      min = 100,
      max = 2000,
      value = 500,
      step = 100,
      ticks = FALSE
    ),

    shiny::plotOutput(
      outputId = "pca_plot",
      height = "600px",
      width = "100%"
    ),


    # Download buttons
    shiny::downloadButton(
      outputId = "download_pca_png",
      label = "Download PNG"
    ),

    shiny::downloadButton(
      outputId = "download_pca_pdf",
      label = "Download PDF"
    )
  ),



  # ===================================
  # TAB 3: VOLCANO PLOT
  # ===================================

  bslib::nav_panel(
    title = "Volcano Plot",

    h2("Differential Gene Expression"),

    p(
      "Explore differentially expressed genes following",
      "SCN5A knockdown in breast cancer samples."
    ),

    # Interactive controls
    bslib::layout_columns(

      shiny::numericInput(
        inputId = "volcano_up_thr",
        label = "Upregulated log2FC threshold",
        value = 1,
        min = 0,
        step = 0.1
      ),

      shiny::numericInput(
        inputId = "volcano_down_thr",
        label = "Downregulated log2FC threshold",
        value = -1,
        max = 0,
        step = 0.1
      ),

      shiny::numericInput(
        inputId = "volcano_p_thr",
        label = "P-value threshold",
        value = 0.05,
        min = 0.0001,
        max = 1,
        step = 0.005
      ),

      shiny::numericInput(
        inputId = "volcano_top_n",
        label = "Top genes to label (per direction)",
        value = 10,
        min = 0,
        max = 30,
        step = 1
      ),

      col_widths = c(3, 3, 3, 3)
    ),

    # Volcano plot
    shiny::plotOutput(
      outputId = "volcano_plot",
      height = "650px",
      width = "100%"
    ),

    # Download buttons
    bslib::layout_columns(

      shiny::downloadButton(
        outputId = "download_volcano_png",
        label = "Download PNG"
      ),

      shiny::downloadButton(
        outputId = "download_volcano_pdf",
        label = "Download PDF"
      ),

      col_widths = c(3, 3)
    )
  ),


  # =====================================
  # Heatmap Analysis Tab
  # =====================================

  bslib::nav_panel(
    title = "Heatmap",

    h2("Differential Gene Expression Heatmap"),

    p(
      "Visualize the expression patterns of the most",
      "upregulated and downregulated genes following",
      "SCN5A knockdown."
    ),

    # Interactive controls
    bslib::layout_columns(

      shiny::numericInput(
        inputId = "heatmap_top_n",
        label = "Top genes per direction",
        value = 10,
        min = 1,
        max = 50,
        step = 1
      ),

      shiny::numericInput(
        inputId = "heatmap_up_thr",
        label = "Upregulated log2FC threshold",
        value = 1,
        step = 0.1
      ),

      shiny::numericInput(
        inputId = "heatmap_down_thr",
        label = "Downregulated log2FC threshold",
        value = -1,
        step = 0.1
      ),

      shiny::numericInput(
        inputId = "heatmap_p_thr",
        label = "P-value threshold",
        value = 0.05,
        min = 0.0001,
        max = 1,
        step = 0.005
      )
    ),

    hr(),

    # Heatmap output
    shiny::plotOutput(
      outputId = "heatmap_plot",
      height = "750px",
      width = "100%"
    ),

    bslib::layout_columns(

      shiny::downloadButton(
        outputId = "download_heatmap_png",
        label = "Download PNG"
      ),

      shiny::downloadButton(
        outputId = "download_heatmap_pdf",
        label = "Download PDF"
      ),

      col_widths = c(3, 3)
    )
  ),


  # =====================================
  # Differential Expression Results Tab
  # =====================================

  bslib::nav_panel(
    title = "DE Results",

    h2("Differential Gene Expression Results"),

    p(
      "Explore DESeq2 results from SCN5A knockdown",
      "versus control breast cancer samples."
    ),

    p(
      "Search for genes, sort columns, and browse",
      "differential expression statistics."
    ),

    hr(),

    # Interactive filtering controls
    bslib::layout_columns(

      shiny::numericInput(
        "de_up_thr",
        "Upregulated log2FC threshold",
        value = 1,
        min = 0,
        step = 0.1
      ),

      shiny::numericInput(
        "de_down_thr",
        "Downregulated log2FC threshold",
        value = -1,
        max = 0,
        step = 0.1
      ),

      shiny::numericInput(
        "de_p_thr",
        "P-value threshold",
        value = 0.05,
        min = 0.0001,
        max = 1,
        step = 0.005
      ),

      shiny::selectInput(
        "de_status",
        "Gene classification",
        choices = c(
          "All",
          "Upregulated",
          "Downregulated",
          "Nonsignificant"
        ),
        selected = "All"
      ),

      col_widths = c(3, 3, 3, 3)
    ),

    hr(),

    h3("Differential Expression Summary"),

    shiny::uiOutput("de_summary_boxes"),

    hr(),

    # DESeq2 results table
    h3("DESeq2 Results Table"),

    # Download button ABOVE the table
    div(
      style = "margin-top: 15px; margin-bottom: 20px;",

      shiny::downloadButton(
        outputId = "download_de_csv",
        label = "Download Filtered Results (CSV)",
        icon = shiny::icon("download"),
        class = "btn-success"
      )
    ),

    # Searchable results table
    DT::DTOutput("de_results_table")

  ),


  # About / Documentation tab
  bslib::nav_panel(
    title = "About",

    h2("About rnaseqViz"),

    p(
      "rnaseqViz is an R Shiny application for exploring",
      "bulk RNA-seq differential gene expression results.",
      "This version uses a breast cancer xenograft dataset",
      "investigating SCN5A knockdown."
    ),

    hr(),

    h3("Dataset Information"),

    tags$ul(
      tags$li(
        strong("GEO Accession: "),
        tags$a(
          "GSE228621",
          href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228621",
          target = "_blank"
        )
      ),
      tags$li(strong("Model: "), "MDA-MB-231 breast cancer xenografts"),
      tags$li(strong("Control: "), "6 samples"),
      tags$li(strong("SCN5A Knockdown: "), "6 samples"),
      tags$li(strong("Genes analyzed: "), "14,703")
    ),

    hr(),

    h3("Differential Expression Analysis"),

    p(
      "Differential gene expression was analyzed using DESeq2.",
      "The application provides gene-level statistics including",
      "log2 fold change, raw p-values, and adjusted p-values."
    ),

    tags$ul(
      tags$li("Default upregulated threshold: log2FC >= 1"),
      tags$li("Default downregulated threshold: log2FC <= -1"),
      tags$li("Default significance threshold: raw p-value < 0.05"),
      tags$li(
        "Thresholds can be adjusted interactively.",
        "Adjusted p-values are provided for interpretation",
        "but are not used for the default classification."
      )
    ),

    hr(),

    h3("Application Features"),

    h4("Principal Component Analysis (PCA)"),

    p(
      "PCA summarizes variation among RNA-seq samples.",
      "Samples are colored by experimental condition",
      "and represented by different point shapes according to cage.",
      "Users can adjust the number of variable genes included."
    ),

    h4("Volcano Plot"),

    p(
      "The volcano plot displays log2 fold change against",
      "-log10(p-value). Genes are classified as upregulated,",
      "downregulated, or nonsignificant based on",
      "user-defined thresholds."
    ),

    h4("Gene Expression Heatmap"),

    p(
      "The heatmap displays the top upregulated and",
      "downregulated genes selected from differential",
      "expression results. Expression values are transformed",
      "and standardized by gene using row-wise Z-scores.",
      "Gene symbols are displayed for interpretation."
    ),

    h4("Differential Expression Results"),

    p(
      "The DE Results tab provides a searchable, sortable,",
      "and filterable table containing gene identifiers,",
      "fold changes, test statistics, and significance values.",
      "Filtered results can be exported as CSV files."
    ),

    hr(),

    h3("Data Visualization and Export"),

    p(
      "PCA plots, volcano plots, and heatmaps can be",
      "downloaded as PNG or PDF files.",
      "Differential expression results can be exported",
      "in CSV format."
    ),

    hr(),

    h3("Software and Technologies"),

    tags$ul(
      tags$li("R"),
      tags$li("Shiny and bslib"),
      tags$li("DESeq2"),
      tags$li("ggplot2"),
      tags$li("pheatmap"),
      tags$li("DT"),
      tags$li("rnaseqViz R package")
    ),

    hr(),

    h3("Source Code"),

    p(
      "The application is developed using reusable",
      "functions implemented in the rnaseqViz R package."
    ),

    p(
      "Source code and installation instructions",
      "will be available through the project's GitHub repository."
    )
  )
)

# Define server

server <- function(input, output, session) {

  # =====================================
  # 1. Reactive PCA plot
  # =====================================

  pca_plot <- shiny::reactive({

    rnaseqViz::plot_pca(
      expression_matrix = expression_matrix,
      metadata = sample_metadata,
      color_by = "condition",
      shape_by = "cage",
      top_n = input$top_n,
      title_text = "PCA of Breast Cancer RNA-seq Samples"
    )

  })

  # =====================================
  # 2. Display PCA plot
  # =====================================

  output$pca_plot <- shiny::renderPlot({

    print(pca_plot())

  }, res = 120)

  # =====================================
  # 3. Download PCA as PNG
  # =====================================

  output$download_pca_png <- shiny::downloadHandler(

    filename = function() {
      paste0(
        "PCA_SCN5A_KD_top",
        input$top_n,
        ".png"
      )
    },

    content = function(file) {

      ggplot2::ggsave(
        filename = file,
        plot = pca_plot(),
        device = "png",
        width = 10,
        height = 7,
        units = "in",
        dpi = 300,
        bg = "white"
      )

    }
  )

  # =====================================
  # 4. Download PCA as PDF
  # =====================================

  output$download_pca_pdf <- shiny::downloadHandler(

    filename = function() {
      paste0(
        "PCA_SCN5A_KD_top",
        input$top_n,
        ".pdf"
      )
    },

    content = function(file) {

      ggplot2::ggsave(
        filename = file,
        plot = pca_plot(),
        device = "pdf",
        width = 10,
        height = 7,
        units = "in",
        bg = "white"
      )

    }
  )



  # =====================================
  # Reactive volcano plot
  # =====================================

  volcano_plot <- shiny::reactive({

    # Validate interactive thresholds
    shiny::req(
      input$volcano_up_thr,
      input$volcano_down_thr,
      input$volcano_p_thr
    )

    shiny::validate(
      shiny::need(
        input$volcano_up_thr > 0,
        "Upregulated threshold must be greater than 0."
      ),

      shiny::need(
        input$volcano_down_thr < 0,
        "Downregulated threshold must be less than 0."
      ),

      shiny::need(
        input$volcano_p_thr > 0 &&
          input$volcano_p_thr <= 1,
        "P-value threshold must be between 0 and 1."
      )
    )

    rnaseqViz::plot_volcano(
      df = res_df,
      title_text = "SCN5A Knockdown vs Control",
      up_thr = input$volcano_up_thr,
      down_thr = input$volcano_down_thr,
      p_thr = input$volcano_p_thr,
      top_n = input$volcano_top_n
    )

  })

  # =====================================
  # Display volcano plot
  # =====================================

  output$volcano_plot <- shiny::renderPlot({

    print(volcano_plot())

  }, res = 120)


  # =====================================
  # Download volcano plot as PNG
  # =====================================

  output$download_volcano_png <- shiny::downloadHandler(

    filename = function() {
      "SCN5A_Volcano_Plot.png"
    },

    content = function(file) {

      ggplot2::ggsave(
        filename = file,
        plot = volcano_plot(),
        device = "png",
        width = 11,
        height = 7,
        units = "in",
        dpi = 300
      )

    }
  )


  # =====================================
  # Download volcano plot as PDF
  # =====================================

  output$download_volcano_pdf <- shiny::downloadHandler(

    filename = function() {
      "SCN5A_Volcano_Plot.pdf"
    },

    content = function(file) {

      ggplot2::ggsave(
        filename = file,
        plot = volcano_plot(),
        device = "pdf",
        width = 11,
        height = 7,
        units = "in"
      )

    }
  )


  # =====================================
  # Heatmap reactive plot
  # =====================================

  heatmap_plot <- shiny::reactive({

    shiny::req(
      input$heatmap_top_n,
      input$heatmap_up_thr,
      input$heatmap_down_thr,
      input$heatmap_p_thr
    )

    shiny::validate(
      shiny::need(
        input$heatmap_up_thr > 0,
        "Upregulated threshold must be positive."
      ),
      shiny::need(
        input$heatmap_down_thr < 0,
        "Downregulated threshold must be negative."
      ),
      shiny::need(
        input$heatmap_p_thr > 0 &&
          input$heatmap_p_thr <= 1,
        "P-value threshold must be between 0 and 1."
      )
    )

    rnaseqViz::plot_heatmap(
      expression_matrix = expression_matrix,
      metadata = sample_metadata,
      de_results = res_df,

      gene_id_col = "entrez_id",
      gene_name_col = "gene_name",

      top_n = input$heatmap_top_n,
      p_thr = input$heatmap_p_thr,
      up_thr = input$heatmap_up_thr,
      down_thr = input$heatmap_down_thr,

      scale_rows = TRUE,
      cluster_rows = TRUE,
      cluster_cols = FALSE,

      annotation_cols = c("condition", "cage"),
      show_rownames = TRUE,

      title_text = "SCN5A Knockdown vs Control"
    )

  })


  # =====================================
  # Render heatmap
  # =====================================

  output$heatmap_plot <- shiny::renderPlot({

    heatmap_plot()

  }, res = 120)


  # =====================================
  # Download Heatmap as PNG
  # =====================================

  output$download_heatmap_png <- shiny::downloadHandler(

    filename = function() {
      "SCN5A_Heatmap.png"
    },

    content = function(file) {

      hm <- heatmap_plot()

      grDevices::png(
        filename = file,
        width = 12,
        height = 9,
        units = "in",
        res = 300
      )

      tryCatch(
        grid::grid.draw(hm$gtable),
        finally = grDevices::dev.off()
      )
    }
  )


  # =====================================
  # Download Heatmap as PDF
  # =====================================

  output$download_heatmap_pdf <- shiny::downloadHandler(

    filename = function() {
      "SCN5A_Heatmap.pdf"
    },

    content = function(file) {

      hm <- heatmap_plot()

      grDevices::pdf(
        file = file,
        width = 12,
        height = 9
      )

      tryCatch(
        grid::grid.draw(hm$gtable),
        finally = grDevices::dev.off()
      )
    }
  )

  # =====================================
  # DE Results: Classify genes
  # =====================================

  de_classified <- shiny::reactive({

    shiny::req(
      input$de_up_thr,
      input$de_down_thr,
      input$de_p_thr
    )

    shiny::validate(
      shiny::need(
        input$de_up_thr > 0,
        "Upregulated threshold must be positive."
      ),
      shiny::need(
        input$de_down_thr < 0,
        "Downregulated threshold must be negative."
      ),
      shiny::need(
        input$de_p_thr > 0 && input$de_p_thr <= 1,
        "P-value threshold must be between 0 and 1."
      )
    )

    df <- as.data.frame(res_df)

    # Classify genes using raw p-values
    up <- !is.na(df$log2FoldChange) &
      !is.na(df$pvalue) &
      df$log2FoldChange >= input$de_up_thr &
      df$pvalue < input$de_p_thr

    down <- !is.na(df$log2FoldChange) &
      !is.na(df$pvalue) &
      df$log2FoldChange <= input$de_down_thr &
      df$pvalue < input$de_p_thr

    df$status <- "Nonsignificant"
    df$status[up] <- "Upregulated"
    df$status[down] <- "Downregulated"

    df
  })


  # =====================================
  # Dynamic DE Summary
  # =====================================

  de_summary <- shiny::reactive({

    df <- de_classified()

    list(
      total = nrow(df),
      up = sum(df$status == "Upregulated"),
      down = sum(df$status == "Downregulated"),
      ns = sum(df$status == "Nonsignificant")
    )
  })


  # =====================================
  # Render Summary Boxes
  # =====================================

  output$de_summary_boxes <- shiny::renderUI({

    stats <- de_summary()

    bslib::layout_columns(

      bslib::value_box(
        title = "Total Genes Analyzed",
        value = format(stats$total, big.mark = ",")
      ),

      bslib::value_box(
        title = "Upregulated",
        value = format(stats$up, big.mark = ",")
      ),

      bslib::value_box(
        title = "Downregulated",
        value = format(stats$down, big.mark = ",")
      ),

      bslib::value_box(
        title = "Nonsignificant",
        value = format(stats$ns, big.mark = ",")
      ),

      col_widths = c(3, 3, 3, 3)
    )
  })


  # =====================================
  # Prepare Filtered DE Table
  # =====================================

  de_table <- shiny::reactive({

    df <- de_classified()

    shiny::req(input$de_status)

    # Filter by gene classification
    if (input$de_status != "All") {
      df <- df[
        df$status == input$de_status,
        ,
        drop = FALSE
      ]
    }

    # Select columns
    df <- df[, c(
      "gene_name",
      "entrez_id",
      "baseMean",
      "log2FoldChange",
      "lfcSE",
      "stat",
      "pvalue",
      "padj",
      "status"
    ), drop = FALSE]

    # Rename columns
    colnames(df) <- c(
      "Gene Symbol",
      "Entrez ID",
      "Base Mean",
      "Log2 Fold Change",
      "LFC Standard Error",
      "Wald Statistic",
      "P-value",
      "Adjusted P-value",
      "Status"
    )

    # Replace missing gene symbols
    missing_symbols <- is.na(df[["Gene Symbol"]]) |
      trimws(as.character(df[["Gene Symbol"]])) == ""

    df[["Gene Symbol"]][missing_symbols] <- "Unknown"

    as.data.frame(df)
  })


  # =====================================
  # Render DE Results Table
  # =====================================

  output$de_results_table <- DT::renderDT({

    df <- de_table()

    DT::datatable(
      df,
      rownames = FALSE,
      filter = "top",
      options = list(
        pageLength = 10,
        lengthMenu = c(10, 25, 50, 100),
        scrollX = TRUE,
        autoWidth = TRUE,
        dom = "lfrtip",
        order = list(list(7, "asc"))
      )
    ) |>
      DT::formatRound(
        columns = c(
          "Base Mean",
          "Log2 Fold Change",
          "LFC Standard Error",
          "Wald Statistic"
        ),
        digits = 3
      ) |>
      DT::formatSignif(
        columns = c(
          "P-value",
          "Adjusted P-value"
        ),
        digits = 3
      )

  }, server = TRUE)


  # Download filtered DESeq2 results
  output$download_de_csv <- shiny::downloadHandler(

    filename = function() {
      paste0(
        "rnaseqViz_DE_results_",
        Sys.Date(),
        ".csv"
      )
    },

    content = function(file) {

      # Retrieve results after classification filtering
      df <- de_table()

      # Retrieve row indices after DT searching
      filtered_rows <- input$de_results_table_rows_all

      if (!is.null(filtered_rows)) {
        df <- df[filtered_rows, , drop = FALSE]
      }

      # Export results
      utils::write.csv(
        x = df,
        file = file,
        row.names = FALSE,
        na = ""
      )
    }
  )

  # Keep CSV download output active
  shiny::outputOptions(
    output,
    "download_de_csv",
    suspendWhenHidden = FALSE
  )


}

# Launch application
shinyApp(ui = ui, server = server)


