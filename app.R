library(shiny)
library(plotly)
library(tibble)
source("./AUCFunction.R")
source("./string_processing.R")

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
      selectInput(
        inputId = 'dataset',
        label = 'Dataset:',
        choices = c(
          'Cortical layers from Developing Human Brain Atlas',
          'Cortical layers from Developing Non-human primate (NHP) Atlas'
        )
      ),
      textAreaInput(
        inputId = "genelist",
        label = "Input your gene list:",
        value = 'MC4R\nADORA1\nZFP179\nGABRB2\nSOX5\nELMO1\nEDIL3\nMGST3',
        rows = 7
      ),
      selectInput(
        inputId = 'species',
        label = 'Species:',
        choices = c('Human', 'Mouse', 'Rhesus Macaque')
      ),
      actionButton(inputId = "submit",
                   label = "Submit"),
      br(),
      br(),
      downloadButton(outputId = "download_data", label = "Download results as .csv")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      div(
        id = "main",
        # Output: Verbatim text for data summary ----
        verbatimTextOutput("summary"),
        br(),
        dataTableOutput("view"),
        br(),
        plotlyOutput("dotplot"),
        verbatimTextOutput("info")
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
    
    cleaned_gene_list <-
      isolate(process_input_genes(input$genelist))
    
    # load reference data
    if (input$dataset == 'Cortical layers from Developing Human Brain Atlas') {
      load('./data/processed/developing_cortical_zones_ranks.Rdata', verbose = TRUE)
      unique_genes <- cortical_zones_ranks$gene_symbol
      paste('loaded data')
      load('./data/processed/developing_cortical_zones_expression_matrix.Rdata',
           verbose = TRUE)
      tidy_expression <- cortical_zones_expression_matrix %>%
        gather(key = zones, value = expression,-gene_symbol)
      
      #target_gene_symbols <- 'Human'
      
      cleaned_gene_list <-
        convert2human(input_genes = cleaned_gene_list, in_species = input$species)
      
    } else {
      load('./data/processed/NHP_cortical_zones_ranks.Rdata', verbose = TRUE)
      cortical_zones_ranks <- NHP_cortical_zones_ranks
      rm(NHP_cortical_zones_ranks)
      unique_genes <- cortical_zones_ranks$gene_symbol
      
      load('./data/processed/NHP_cortical_zones_expression_matrix.Rdata',
           verbose = TRUE)
      tidy_expression <- NHP_cortical_zones_expression_matrix %>%
        gather(key = zones, value = expression,-gene_symbol)
      
      #target_gene_symbols <- 'NHP'
      
      cleaned_gene_list <-
        convert2nhp(input_genes = cleaned_gene_list, in_species = input$species)
    }
    
    #cleaned_gene_list <- isolate(process_input_genes(input$genelist))
    #if (input$species == 'Mouse') {
    #  cleaned_gene_list <- convert_genes(cleaned_gene_list)
    #}
    # print to console
    print(paste0("Before time taken:", Sys.time() - start))
    
    #for indices - use dplyr for ease
    forIndices <- as_tibble(cortical_zones_ranks$gene_symbol)
    names(forIndices) <- 'gene_symbol'
    forIndices %<>% mutate(isTargetGene = gene_symbol %in% cleaned_gene_list)
    targetIndices <- forIndices$isTargetGene
    
    # only columns from cortical zones remain in df
    df <- cortical_zones_ranks %>%
      select(-gene_symbol)
    
    #AUROC <- map_df(df, auroc_analytic, as.numeric(targetIndices))
    AUROC <- map_df(df, auroc_analytic, targetIndices)
    wilcox_tests <- map_df(df, apply_MWU, targetIndices)
    
    # group results together in a single table
    table <-
      bind_cols(gather(AUROC, key = zone, value = AUROC),
                gather(wilcox_tests, value = pValue)) %>%
      select(-key)
    
    print(paste0("Wilcox time taken:", Sys.time() - start))
    
    # these are the values for the results table
    table %<>% arrange(-AUROC)
    table %<>% mutate(
      pValue = signif(pValue, digits = 3),
      AUROC = signif(AUROC, digits = 3),
      adjusted_P = signif(p.adjust(pValue), digits = 3)
    )
    #table$order
    x <- tibble(order = 1:7, zone = c("ventricular zone", "subventricular zone", "intermediate zone", "subplate zone", "cortical plate",
                                      "marginal zone","subpial granular zone"))
    table <- inner_join(table, x, by = 'zone')
    
    selected_values <- reactive({
      req(cleaned_gene_list)
      selected_values <-
        tidy_expression %>% filter(gene_symbol %in% cleaned_gene_list)
      if (input$dataset == 'Cortical layers from Developing Human Brain Atlas') {
        selected_values %<>% mutate(zones = factor(zones, levels = c("ventricular zone", "subventricular zone", "intermediate zone", "subplate zone",
                                                                     "cortical plate", "marginal zone", "subpial granular zone")))
      } else {
        # input$dataset == 'Cortical layers from Developing Non-human primate (NHP) Atlas'
        selected_values %<>% mutate(zones = factor(zones, levels = c("ventricular zone", "ventricular zone (inner)", "ventricular zone (outer)",
                                                                     "subventricular zone", "subventricular zone (inner)", "subventricular zone (outer)",
                                                                     "intermediate zone", "inner fiber zone", "outer fiber zone", "transitory migratory zone",
                                                                     "subplate zone", "cortical plate", "cortical plate (inner)", "cortical plate (outer)",
                                                                     "marginal zone", "subpial granular zone", "white matter", "layer I", "layer II", 
                                                                     "layer II/III", "layer III", "layer IV", "layer V", "layer VI")))
      }
    })
    
    # these are the values used for the plot
    #selected_values <- tidy_expression %>% filter(gene_symbol %in% cleaned_gene_list)
    #selected_values %<>% mutate(zones = factor(zones, levels=c("ventricular zone", "subventricular zone", "intermediate zone", "subplate zone", "cortical plate", "marginal zone", "subpial granular zone")))
    
    output$summary <- renderPrint({
      #count of intersection of submitted genes with total gene list
      cat(paste("Time taken:", round(Sys.time() - start), "seconds"))
      cat(paste(
        "\nGenes found in data:",
        sum(cleaned_gene_list %in% unique_genes),
        "of",
        length(cleaned_gene_list)
      ))
    })
    
    output$view <- renderDataTable({
      table
    }, escape = FALSE)
    
    output$dotplot <- renderPlotly({
      if (length(cleaned_gene_list) > 20) {
        p <- ggplot(selected_values(), aes(x = zones, y = expression)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(alpha = 0.25, width = 0.1) +
          theme(axis.text.x = element_text(angle = 45, hjust = 0))
        
        ggplotly(p) #%>% layout(dragmode = "select")
        
      } else {
        p <- ggplot(selected_values(), aes(x = zones, y = expression)) +
          #geom_dotplot(binaxis = "y", stackdir = "center")
          geom_jitter(width = 0.1) +
          theme(axis.text.x = element_text(angle = 45, hjust = 0))
        
        ggplotly(p) #%>% layout(dragmode = "select")
      }
    })
    
    output$download_data <-
      downloadHandler(
        filename = "polygenic_layers_AUC_results.csv",
        content = function(file) {
          write_csv(table, file)
        }
      )
    
  })
}

shinyApp(ui, server)
