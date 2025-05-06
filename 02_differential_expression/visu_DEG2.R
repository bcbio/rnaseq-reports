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
import::from(data.table,fread)
inputRead <- function(f){
  if(sum(endsWith(f,c("rds","RDS")))>0){
    if(R.utils::isUrl(f)){f <- url(f)}
    return(readRDS(f))
  }else if(sum(endsWith(f,c("qs","QS")))>0){
    if(R.utils::isUrl(f)){f <- url(f)}
    return(qread(f))
  }else{
    print("Check file extension and choose appropriate functions!")
  }
}

# the folder containing all your DEG results
DEG_outputF <- paste0("https://raw.githubusercontent.com/",
                      "bcbio/bcbioR-test-data/main/rnaseq/DEG_visualization")
      
# the factor of interest in your data that you are trying to evaluate
column <- "sample_type"
# contrast groups
contrasts <- "normal_vs_tumor"

deg_pathway <- glue("{DEG_outputF}/full_{column}_{contrasts}_pathways.csv")
deg_result <-  glue("{DEG_outputF}/full_{column}_{contrasts}_deg.csv")
full_expr <- glue("{DEG_outputF}/full_expression.csv")
metaD <- glue("{DEG_outputF}/deseq_coldata.rds")


deg_annotated <- fread(deg_pathway) %>% filter(padj < 0.05)
gene_list <- map(deg_annotated$genes,\(d) 
                 setdiff(strsplit(d,split=",")[[1]],"")) %>% 
  unlist() %>% unique()
volcano_data <- fread(deg_result) %>%
  select(gene = gene_name, logFC = lfc, adj.P.Val = padj) %>% 
  filter(gene!="")

expression_data <- fread(full_expr) %>%
  select(gene = gene_name, starts_with("SRX")) %>%
  filter(gene %in% gene_list)

meta <- inputRead(metaD) %>%
  select(sample, condition = sample_type) %>%
  magrittr::set_rownames(NULL) %>%
  arrange(condition)

ui <- dashboardPage(
  dashboardHeader(title = "Gene-Centric DEG View"),
  dashboardSidebar(
    selectInput("selected_gene", "Select Gene", choices = sort(gene_list))
  ),
  dashboardBody(
    tags$head(tags$style(HTML(".tab-content { padding-top: 10px; }"))),
    fluidRow(
      box(title = "Pathways Associated with Gene", width = 12, solidHeader = TRUE, status = "primary",
          DTOutput("pathway_table"))
    ),
    fluidRow(
      column(width = 6,
             box(title = "Volcano Plot", width = NULL, solidHeader = TRUE, status = "primary",
                 plotlyOutput("volcano_plot"))
      ),
      column(width = 6,
             box(title = "Expression Boxplot", width = NULL, solidHeader = TRUE, status = "primary",
                 plotOutput("expression_plot"))
      )
    ),
    tags$style(HTML("
      .content-wrapper, .right-side {
        overflow-x: auto;
      }
      .content {
        min-width: 1200px;
      }
    "))
  ),
  skin = "blue"
)

server <- function(input, output, session) {
  
  output$pathway_table <- renderDT({
    req(input$selected_gene)
    g <- input$selected_gene
    genePat <- glue("(?<=^|,){g}(?=,|$)")
    
    deg_annotated %>%
      filter(grepl(genePat, genes, perl = TRUE)) %>%
      select(-genes, -comparison) %>%
      arrange(padj) %>% 
      mutate(
        padj = formatC(padj, format = "E", digits = 2),
        NES = formatC(NES, digits = 3)
      ) %>% 
      datatable(options = list(pageLength = 10))
  })
  
  output$volcano_plot <- renderPlotly({
    req(input$selected_gene)
    df <- volcano_data
    df$highlight <- ifelse(df$gene == input$selected_gene, "Selected", "Other")
    
    p <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val),
                        text = paste("Gene:", gene,
                                     "<br>logFC:", round(logFC, 2),
                                     "<br>adj.P.Val:", signif(adj.P.Val, 3)))) +
      geom_point(aes(color = highlight)) +
      scale_color_manual(values = c("Selected" = "darkred", "Other" = "gray")) +
      theme_minimal()
    
    ggplotly(p, tooltip = "text")
  })
  
  output$expression_plot <- renderPlot({
    req(input$selected_gene)
    g <- input$selected_gene
    
    barplotD <- subset(expression_data, gene == g) %>%
      tidyr::gather(sample, expression, -gene) %>%
      inner_join(meta, by = "sample")
    
    ggpubr::ggboxplot(barplotD,
                      x = "condition", y = "expression",
                      color = "condition", palette = "jco", 
                      title = g,
                      add = "jitter",
                      add.params = list(size = rel(3), alpha = 0.6)) +
      theme_minimal()+
      theme(axis.text.x = element_text(size = rel(2)),
            axis.title.y = element_text(size = rel(2)),
            plot.title = element_text(size = rel(2.2), face = "bold", hjust = 0.5),
            legend.position="none")+
      labs(x="")
  })
}

shinyApp(ui, server)