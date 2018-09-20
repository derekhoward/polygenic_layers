library(shiny)

source("./AUCFunction.R")
source("./string_processing.R")

load('./processed_developmental_zones.Rdata', verbose=TRUE)
unique_genes <- result$gene_symbol

apply_MWU <- function(column, targetIndices) {
  wilcox.test(column[targetIndices], column[!targetIndices], conf.int = F)$p.value
}

ui <- fluidPage(
  # App title ----
  titlePanel("Polygenic tester for the developing human brain"),
  
  # Sidebar layout with input and output definitions ----
    # Sidebar panel for inputs ----
    sidebarPanel(
      # Input: Selector for choosing dataset ----
      textAreaInput(inputId = "genelist",
                    label = "Input your gene list:",
                    value = 'MC4R\nADORA1\nZFP179\nGABRB2\nSOX5\nELMO1\nEDIL3\nMGST3',
                    rows=7),
      selectInput('species', 'Species:',
                  choices=c('Human', 'Mouse')),
      actionButton("submit", "Submit")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      div(id = "main",
          
      # Output: Verbatim text for data summary ----
      verbatimTextOutput("summary"),
      br(),
      #plotOutput("plot"),
      br(),
      # Output: HTML table with requested number of observations ----
      dataTableOutput("view")
      
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
    forIndices <- as_tibble(result$gene_symbol)
    names(forIndices) <- 'gene_symbol'
    forIndices %<>% mutate(isTargetGene = gene_symbol %in% cleaned_gene_list)
    targetIndices <- forIndices$isTargetGene
    
    # only columns from cortical zones remain in df
    df <- result %>% 
      select(-gene_symbol)
    
    AUROC <- map_df(df, auroc_analytic, as.numeric(targetIndices))
    wilcox_tests <- map_df(df, apply_MWU, targetIndices)
    
    # group results together in a single table
    table <- bind_cols(gather(AUROC, key = zone, value = AUROC), gather(wilcox_tests, value = pValue)) %>% 
      select(-key)
    
    print(paste0("Wilcox time taken:", Sys.time() - start))
    
    #add descriptions
    table %<>% arrange(-AUROC)
    
    output$summary <- renderPrint({
      #count of intersection of submitted genes with total gene list
      cat(paste("Time taken:", round(Sys.time() - start), "seconds"))
      cat(paste("\nGenes found in data:",sum(cleaned_gene_list %in% unique_genes), "of", length(cleaned_gene_list)))
    })
    
    output$view <- renderDataTable({
      table %<>% mutate(pValue = signif(pValue, digits=3), AUROC = signif(AUROC, digits=3), adjusted_P = signif(p.adjust(pValue), digits=3))
      table
    }, escape = FALSE)
  }
  )
}

shinyApp(ui, server)

