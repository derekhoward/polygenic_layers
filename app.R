library(shiny)
library(ggplot2)
library(magrittr)
library(plotly)
library(tibble)
library(tidyr)
library(dplyr)
library(purrr)
library(shinyjs)
source("./AUCFunction.R")
source("./string_processing.R")

apply_MWU <- function(column, targetIndices) {
  wilcox.test(column[targetIndices], column[!targetIndices], conf.int = F)$p.value
}

ui <- fluidPage(
  shinyjs::useShinyjs(),
  tags$head(includeHTML("google-analytics.html")),
  # App title ----
  titlePanel("Polygenic tester for developing cortical layers"),
  
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
        value = 'THRA\nRTN1\nTUBA1A\nSTMN2\nCRMP1\nTUBB3\nISLR2',
        rows = 10
      ),
      selectInput(
        inputId = 'species',
        label = 'Species of input genes:',
        choices = c('Human', 'Mouse', 'Rhesus Macaque')
      ),
      actionButton(inputId = "submit",
                   label = "Submit"),
      br(),
      br(),
      downloadButton(outputId = "download_data", label = "Download results as .csv"),
      hr(),
      tags$b("Data was made available by the Allen Institute for Brain Science and is available from: "),
      br(),
      tags$a(href="http://www.brainspan.org/static/download.html", "Developing Human Brain Atlas"),
      br(),
      tags$a(href="http://www.blueprintnhpatlas.org/static/download", "Cortical layers from Developing Non-human primate (NHP) Atlas")
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
    cat("\nResults will load here when complete")
    cat("\n")
    #print(gc())
    #print(Sys.info()['nodename'])
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
      load('./data/processed/developing_cortical_zones_expression_matrix.Rdata', verbose = TRUE)
      tidy_expression <- cortical_zones_expression_matrix %>%
        gather(key = zones, value = expression,-gene_symbol)
      
      #target_gene_symbols <- 'Human'
      
      cleaned_gene_list <- convert2human(input_genes = cleaned_gene_list, in_species = input$species)
      
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
    
    AUROC <- map_df(df, auroc_analytic, targetIndices)
    wilcox_tests <- map_df(df, apply_MWU, targetIndices)
    
    # group results together in a single table
    table <- bind_cols(gather(AUROC, key = zone, value = AUROC), 
                       gather(wilcox_tests, value = pValue)) %>%
      select(-key)
    
    print(paste0("Wilcox time taken:", Sys.time() - start))
    
    # these are the values for the results table
    table %<>% arrange(-AUROC)
    table %<>% mutate(pValue = signif(pValue, digits = 3), 
                      AUROC = signif(AUROC, digits = 3),
                      adjusted_P = signif(p.adjust(pValue), digits = 3))
    
    #table$order
    if (input$dataset == 'Cortical layers from Developing Human Brain Atlas') {
      zone_order <-c("ventricular zone", "subventricular zone", "intermediate zone", "subplate zone", "cortical plate",
                     "marginal zone","subpial granular zone")
      x <- tibble(laminar_order = 1:7, zone = zone_order)
    } else {
      #from Blueprint documentation: http://download.alleninstitute.org/nhp/Prenatal_Macaque_LMD_Microarray/neuroanatomical_guides_for_LMD_sampling/V1/
      #The layers annotated for LMD are indicated in the right panel and include marginal zone (mz), layer 2 (2), layer 2/3 (2/3), layer 3 (3), 
      zone_order <- c("white matter", "marginal zone", "layer I", "layer II", "layer II/III", "layer III", 
                      #outer cortical plate (cpo), layer 4 (4), layer 4A (4A), layer 4B (4B), layer 4Ca (4Ca), layer 4Cb (4Cb), layer 5 (5), layer 6 (6), cortical plate (cp), 
                      "cortical plate (outer)",  "layer IV", "layer V", "layer VI", "cortical plate", 
                      #inner cortical plate (cpi), subplate (sp), intermediate zone (iz), intermediate cell dense zone (icd), transitory migratory zone (tmz), 
                      "cortical plate (inner)","subplate","intermediate zone","transitory migratory zone", 
                      #outer fiber (plexiform) zone (ofz), subventricular zone (sz), outer subventricular zone (szo), inner fiber (plexiform) zone (ifz), 
                      "outer fiber zone","subventricular zone","subventricular zone (outer)", "inner fiber zone", 
                      #inner subventricular zone (szi), outer ventricular zone (vzo), inner ventricular zone (vzi), and ventricular zone (vz). 
                      "subventricular zone (inner)","ventricular zone (outer)","ventricular zone (inner)", "ventricular zone")
      zone_order <- rev(zone_order) #reverse so it lines up with human direction
      
      x <- tibble(laminar_order = 1:length(zone_order), zone = zone_order)
   }
    table <- inner_join(table, x, by = 'zone')

    selected_values <- reactive({
      req(cleaned_gene_list)
      selected_values <- tidy_expression %>% filter(gene_symbol %in% cleaned_gene_list)
      selected_values %<>% mutate(zones = factor(zones, levels = zone_order))
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
      #base plot
      p <- ggplot(selected_values(), aes(x = zones, y = expression, names=gene_symbol))
      #figure depends on size
      if (length(cleaned_gene_list) > 20) {
        p <-  p + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.25, width = 0.1)
      } else {
        p <- p + geom_jitter(width = 0.1)
      }
      p <- p + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 0))
      p <- p + geom_hline(yintercept = 0, color='darkgrey', size=0.4)
      ggplotly(p) #%>% layout(dragmode = "select")
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
