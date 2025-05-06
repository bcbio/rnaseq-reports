library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(ComplexHeatmap)
library(ggpubr)
library(grid)
library(purrr)
library(dplyr)
library(glue)
import::from(data.table, fread)
inputRead <- function(f) {
  if (sum(endsWith(f, c("rds", "RDS"))) > 0) {
    if (R.utils::isUrl(f)) {
      f <- url(f)
    }
    return(readRDS(f))
  } else if (sum(endsWith(f, c("qs", "QS"))) > 0) {
    if (R.utils::isUrl(f)) {
      f <- url(f)
    }
    return(qread(f))
  } else {
    print("Check file extension and choose appropriate functions!")
  }
}

# the folder containing all your DEG results
DEG_outputF <- "https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/DEG_visualization"
# the factor of interest in your data that you are trying to evaluate
column <- "sample_type"
# contrast groups
contrasts <- "normal_vs_tumor"

deg_pathway <- glue("{DEG_outputF}/full_{column}_{contrasts}_pathways.csv")
deg_result <- glue("{DEG_outputF}/full_{column}_{contrasts}_deg.csv")
full_expr <- glue("{DEG_outputF}/full_expression.csv")
metaD <- glue("{DEG_outputF}/deseq_coldata.rds")


deg_annotated <- fread(deg_pathway) %>% filter(padj < 0.05)
deg_lists <- split(deg_annotated$genes, deg_annotated$pathway)
deg_lists <- map(deg_lists, \(d) setdiff(strsplit(d, split = ",")[[1]], ""))
geneset_list <- deg_lists[lengths(deg_lists) >= 5]

volcano_data <- fread(deg_result) %>%
  select(gene = gene_name, logFC = lfc, adj.P.Val = padj)

expression_data <- fread(full_expr) %>%
  select(gene = gene_name, starts_with("SRX")) %>%
  filter(gene %in% unique(unlist(geneset_list)))
meta <- inputRead(metaD) %>%
  select(sample, condition = sample_type) %>%
  magrittr::set_rownames(NULL) %>%
  arrange(condition)
heatmap_matrix <- as.matrix(expression_data %>% tibble::column_to_rownames("gene"))[, meta$sample]
colPal <- c("#1f77b4", "#ff7f0e")
conditions <- levels(meta$condition)
col_annot <- HeatmapAnnotation(
  Condition = setNames(meta$condition, meta$sample),
  col = list(Condition = setNames(colPal[1:length(conditions)], conditions)),
  show_annotation_name = TRUE,
  annotation_legend_param = list(title = "", ncol = 2)
)


ui <- dashboardPage(
  dashboardHeader(title = "DEG Pathway"),
  dashboardSidebar(
    selectInput("geneset", "Select Gene Set", choices = names(geneset_list))
  ),
  dashboardBody(
    tags$head(tags$style(HTML(".tab-content { padding-top: 10px; }"))),
    fluidRow(
      column(
        width = 4,
        box(
          title = "Genes in Gene Set", width = NULL, solidHeader = TRUE, status = "primary",
          DTOutput("gene_table")
        )
      ),
      column(
        width = 4,
        box(
          title = "Volcano Plot", width = NULL, solidHeader = TRUE, status = "primary",
          plotlyOutput("volcano_plot")
        ),
        box(
          title = "Expression Boxplots", width = NULL, solidHeader = TRUE, status = "primary",
          # Replace tabs with a dropdown selector for genes
          selectInput("selected_gene", "Select Gene", choices = NULL),
          plotOutput("expression_plot")
        )
      ),
      column(
        width = 4,
        box(
          title = "Heatmap of Selected Genes", width = NULL, solidHeader = TRUE, status = "primary",
          plotOutput("heatmap")
        )
      ),
      tags$style(HTML("
        .content-wrapper, .right-side {
          overflow-x: auto;
        }
        .content {
          min-width:1500px;
        }
      "))
    )
  ),
  skin = "blue" # Optional theme color
)

# Server
server <- function(input, output, session) {
  selected_genes <- reactive({
    geneset_list[[input$geneset]]
  })

  # Update gene dropdown when geneset changes
  observe({
    genes <- selected_genes()
    if (length(genes) > 0) {
      updateSelectInput(session, "selected_gene", choices = genes, selected = genes[1])
    } else {
      updateSelectInput(session, "selected_gene", choices = character(0))
    }
  })

  output$gene_table <- renderDT({
    genetab <- volcano_data %>%
      filter(gene %in% selected_genes()) %>%
      mutate(adj.P.Val = formatC(adj.P.Val, format = "E", digits = 3)) %>%
      mutate(logFC = formatC(logFC, digits = 3))
    datatable(genetab, selection = "multiple")
  })

  output$volcano_plot <- renderPlotly({
    df <- volcano_data
    df$highlight <- ifelse(df$gene %in% selected_genes(), "Selected", "Other")

    p <- ggplot(df, aes(
      x = logFC, y = -log10(adj.P.Val),
      text = paste(
        "Gene:", gene,
        "<br>logFC:", round(logFC, 2),
        "<br>adj.P.Val:", round(adj.P.Val, 4)
      )
    )) +
      geom_point(aes(color = highlight)) +
      scale_color_manual(values = c("Selected" = "darkred", "Other" = "gray")) +
      theme_minimal()

    ggplotly(p, tooltip = "text")
  })

  # Render boxplot for the selected gene
  output$expression_plot <- renderPlot({
    req(input$selected_gene)
    g <- input$selected_gene
    barplotD <- subset(expression_data, gene == g) %>%
      tidyr::gather(sample, expression, -gene) %>%
      inner_join(meta, by = "sample")

    ggpubr::ggboxplot(barplotD,
      x = "condition", y = "expression",
      color = "condition", palette = "jco", title = g
    ) +
      theme_minimal()
  })

  output$heatmap <- renderPlot({
    genes <- selected_genes()
    gene_indices <- rownames(heatmap_matrix) %in% genes
    mat <- heatmap_matrix[gene_indices, , drop = FALSE]

    inferno_colors <- viridis::inferno(10)

    ht <- Heatmap(mat,
      name = "Expression",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      # Add these parameters for customization:
      row_names_side = "left", # Display row names on the left
      column_names_rot = 45, # Rotate column names by 45 degrees
      col = inferno_colors, # Use inferno color palette with 10 colors
      heatmap_legend_param = list(title = "Expression", ncol = 1),
      top_annotation = col_annot
    )


    grid.newpage()
    draw(ht, heatmap_legend_side = "right", annotation_legend_side = "top")
  })
}

shinyApp(ui, server)
