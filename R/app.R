library(shiny)
library(ggplot2)
source("./AUCFunction.R")
source("./string_processing.R")

#loads ranked genes in different layers as 'cortical_zones_ranks'
load('../data/developmental_zones_ranks.Rdata', verbose=TRUE)
unique_genes <- cortical_zones_ranks$gene_symbol

#loads zscored genes in different layers as 'cortical_zones_expression_matrix'
load('../data/cortical_zones_expression_matrix.Rdata', verbose=TRUE)
tidy_expression <- cortical_zones_expression_matrix %>% 
  gather(key = zones, value = expression, -gene_symbol)

apply_MWU <- function(column, targetIndices) {
  wilcox.test(column[targetIndices], column[!targetIndices], conf.int = F)$p.value
}

ui <- fluidPage(
  # App title ----
  titlePanel("Polygenic tester for the developing human brain"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      # Input: Selector for choosing dataset ----
      textAreaInput(inputId = "genelist",
                    label = "Input your gene list:",
                    value = 'MC4R\nADORA1\nZFP179\nGABRB2\nSOX5\nELMO1\nEDIL3\nMGST3',
                    rows=7),
      selectInput(inputId = 'species', 
                  label = 'Species:',
                  choices=c('Human', 'Mouse')),
      actionButton(inputId = "submit", 
                   label = "Submit")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      div(id = "main",
          # Output: Verbatim text for data summary ----
          verbatimTextOutput("summary"),
          br(),
          dataTableOutput("view"),
          br(),
          #plotOutput("dotplot", hover='plot_hover'),
          plotOutput("dotplot", brush='plot_brush'),
          verbatimTextOutput("info"),
          downloadButton(outputId = "download_data", label = "Download data")
    )
  )
  )
)


# Define server logic process and output top cortical layers/zones ----
server <- function(input, output) {
  output$summary <- renderPrint({
    #cat(paste("cores set to", cores))
    cat("\nResults will load here when complete")
    cat("\n")
    print(gc())
    print(Sys.info()['nodename'])
  })
  
  observeEvent(input$submit, {
    start <- Sys.time()
    cleaned_gene_list <- isolate(process_input_genes(input$genelist))
    if (input$species == 'Mouse') {
      cleaned_gene_list <- convert_genes(cleaned_gene_list)
    }
    print(paste0("Before time taken:", Sys.time() - start))
    
    #for indices - use dplyr for ease
    forIndices <- as_tibble(cortical_zones_ranks$gene_symbol)
    names(forIndices) <- 'gene_symbol'
    forIndices %<>% mutate(isTargetGene = gene_symbol %in% cleaned_gene_list)
    targetIndices <- forIndices$isTargetGene
    
    # only columns from cortical zones remain in df
    df <- cortical_zones_ranks %>% 
      select(-gene_symbol)
    
    AUROC <- map_df(df, auroc_analytic, as.numeric(targetIndices))
    wilcox_tests <- map_df(df, apply_MWU, targetIndices)
    
    # group results together in a single table
    table <- bind_cols(gather(AUROC, key = zone, value = AUROC), gather(wilcox_tests, value = pValue)) %>% 
      select(-key)
    
    print(paste0("Wilcox time taken:", Sys.time() - start))
    
    #add descriptions
    table %<>% arrange(-AUROC)
    
    selected_values <- tidy_expression %>% filter(gene_symbol %in% cleaned_gene_list)
    
    output$summary <- renderPrint({
      #count of intersection of submitted genes with total gene list
      cat(paste("Time taken:", round(Sys.time() - start), "seconds"))
      cat(paste("\nGenes found in data:",sum(cleaned_gene_list %in% unique_genes), "of", length(cleaned_gene_list)))
    })
    
    output$view <- renderDataTable({
      table %<>% mutate(pValue = signif(pValue, digits=3), AUROC = signif(AUROC, digits=3), adjusted_P = signif(p.adjust(pValue), digits=3))
      table
    }, escape = FALSE)
    
    output$dotplot <- renderPlot({
      if (length(cleaned_gene_list) > 20) {
        ggplot(selected_values, aes(x=zones, y=expression)) +
          geom_boxplot()
      } else {
        ggplot(selected_values, aes(x=zones, y=expression)) +
          geom_dotplot(binaxis = "y", stackdir = "center")
      }
    })
    
    output$info <- renderPrint({
      #nearPoints(selected_values, coordinfo = input$plot_hover, xvar = "zones", yvar = "expression")
      brushedPoints(selected_values, brush = input$plot_brush, xvar = "zones", yvar = "expression")
    })
    
    output$download_data <- downloadHandler(filename = "polygenic_layers_AUC_results.csv", 
                                            content = function(file) { write_csv(table, file) })
      
  }
  )
}

shinyApp(ui, server)

