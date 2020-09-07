library(shiny)
library(tidyverse)
library(shinydashboard)
library(data.table)
library(tximport)
library(DESeq2)
library(IHW)
library(plotly)
library(cowplot)
library(ggsci)

meta_trap <- readRDS("data/meta_trap.rds") %>%
  filter(ip_pool == "n")
setDT(meta_trap)
tx2gene <- readRDS("data/tx2gene_trap.rds")

#Custom functions
#Make tidy data.frame
make_tidy <- function(x, first_column){
  as.data.frame(x) %>%
    rownames_to_column(first_column)
}
#summarise res
make_res_summary <- function (res, fdr) {
  res %>%
    make_tidy("gene") %>%
    filter(padj < fdr) %>%
    mutate(Sign = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated")) %>%
    group_by(Sign) %>%
    summarise(Number = n())
}

header <- dashboardHeader(title = "isitdifferentiallyexpressed?")

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Step 1: Choose the Dataset", tabName = "metadata_build", icon = icon("table")), 
    menuItem("Step 2: Produce the Results", tabName = "results", icon = icon("check-circle"))
  )
)
 
 body <- dashboardBody(
  tabItems(
    tabItem(tabName = "metadata_build",
            fluidRow(
              column(width = 3, 
                     box(title = "Introduction", 
                         height = NULL, 
                         collapsible = FALSE, 
                         solidHeader = TRUE, 
                         status = "info", 
                         width = NULL, 
                         p("This is an introduction to the website.")
                         ), 
                     selectInput(inputId = "experiment_choice", 
                                 label = "Choose Experiment", 
                                 choices = c("TRAP", "iPSC"), 
                                 selected = NULL,
                                 selectize = TRUE), 
                     conditionalPanel(
                       condition = 'input.experiment_choice == "TRAP"', 
                       selectizeInput(inputId = "trap_choice", 
                                     label = "Which TRAP Cohort?", 
                                     choices = unique(meta_trap$cohort))),
                     conditionalPanel(
                       condition = 'input.experiment_choice == "TRAP"', 
                       checkboxGroupInput(inputId = "age_choice", 
                                      label = "Which age groups?", 
                                      choices = unique(meta_trap$age), 
                                      select = unique(meta_trap$age))),
                     conditionalPanel(
                       condition = 'input.trap_choice != "PK1"', 
                       checkboxGroupInput(inputId = "ovx_choice", 
                                      label = "Which genotypes?", 
                                      choices = unique(meta_trap$ovx), 
                                      select = unique(meta_trap$ovx))),
                     conditionalPanel(
                       condition = 'input.experiment_choice == "TRAP" & input.trap_choice == "PK1"', 
                       checkboxGroupInput(inputId = "region_choice", 
                                      label = "Which regions?", 
                                      choices = unique(meta_trap$region), 
                                      select = unique(meta_trap$region))),
                     conditionalPanel(
                       condition = 'input.experiment_choice == "TRAP"', 
                       checkboxInput(inputId = "ip_only",
                                     label = "TRAP Samples Only", 
                                     value = TRUE)),
                     conditionalPanel(
                       condition = 'input.trap_choice == "PK1"', 
                       checkboxInput(inputId = "remove_outliers", 
                                     label = "Remove outliers", 
                                     value = TRUE)),
                     selectizeInput(inputId = "exclude_samples_choice", 
                                    label = "Exclude any samples?", 
                                    choices = NULL, 
                                    multiple = TRUE), 
                     selectizeInput(inputId = "comparison_choice", 
                                    label = "Comparison Factor", 
                                    choices = NULL),
                     fluidRow(
                       column(
                         width = 4, 
                         selectizeInput(inputId = "nominator_choice", 
                                        label = "Nominator", 
                                        choices = NULL)
                       ), 
                       column(
                         width = 1, 
                         hr(),
                         tags$b("vs")
                       ),
                       column(
                         width = 4, 
                         selectizeInput(inputId = "denominator_choice", 
                                        label = "Denominator", 
                                        choices = NULL)
                       )
                     ),
                     selectizeInput(inputId = "covariates_choice", 
                                    label = "Covariates to include:", 
                                    multiple = TRUE,
                                    choices = NULL),
                     verbatimTextOutput("final_rows_value"),
                     uiOutput("ui_makedds"), 
                     hr(),
                     verbatimTextOutput("debugdiy"), 
                     box(
                       title = "Result Summary", 
                       DT::dataTableOutput("res_summary"), 
                       width = NULL, 
                       solidHeader = TRUE
                     )
                     ),
              column(
                width = 9, 
                box(
                  title = "TRAP Metadata",
                  DT::dataTableOutput("meta_select"), 
                  width = NULL, 
                  solidHeader = TRUE
                  ), 
                box(title = "Results Plot", 
                    plotlyOutput("results_plot"))
                )
              )
            ),
    
    
    #NEXT PAGE
    tabItem(
      tabName = "results", 
      fluidRow(
        column(width = 3, 
               box(title = "Introduction", 
                   height = NULL, 
                   collapsible = FALSE, 
                   solidHeader = TRUE, 
                   status = "info", 
                   width = NULL, 
                   p("This is the results page.")
               )
      )
      )
    )
    )
  )

ui <- dashboardPage(header, 
                    sidebar, 
                    body)

server <- function(input, output, session) {
  
  #Make a list of reactive values - not currently in use
  values <- reactiveValues()
  
  #Rows selected by TRAP choice
  trap_choice_rows <- reactive({
    meta_trap[cohort %in% input$trap_choice, which = TRUE]
  })
  #Rows selected by excluding outliers
  outlier_rows <- reactive({
    if (input$remove_outliers == TRUE) {
    meta_trap[outlier == FALSE, which = TRUE]
    } else {
    meta_trap[outlier == TRUE | outlier == FALSE, which = TRUE]
    }
  })
  #Rows selected by choosing only IP samples
  ip_rows <- reactive({
    if (input$ip_only == TRUE) {
      meta_trap[ip == "ip", which = TRUE]
    } else {
      meta_trap[ip %in% meta_trap$ip, which = TRUE]
    }
  })
  #Rows from age_choice
  age_choice_rows <- reactive({
    meta_trap[age %in% input$age_choice, which = TRUE]
  })
  #Rows from OVX choice
  ovx_choice_rows <- reactive({
    meta_trap[ovx %in% input$ovx_choice, which = TRUE]
  })
  #Rows from region choice
  region_choice_rows <- reactive({
    meta_trap[region %in% input$region_choice, which = TRUE]
  })
  
  meta0 <- reactive({
    rows <- Reduce(intersect, 
                   list(trap_choice_rows(),
                        outlier_rows(),
                        ip_rows(), 
                        age_choice_rows(), 
                        ovx_choice_rows(), 
                        region_choice_rows()))
    meta_trap[rows]
  })
  observeEvent(meta0(), {
    updateSelectizeInput(session,
                         "exclude_samples_choice", 
                         choices = meta0()$name, 
                         selected = lapply(reactiveValuesToList(input), unclass)$exclude_samples_choice)
  })
  
  #Rows excluded specifically
  exclude_samples_rows <- reactive({
    meta_trap[!name %in% input$exclude_samples_choice, which = TRUE ]
  })
  
  
  #Build the filtered table by intersecting a list of all the row numbers
  meta <- reactive({
    final_rows <- Reduce(intersect,
                         list(trap_choice_rows(),
                              outlier_rows(),
                              ip_rows(), 
                              age_choice_rows(), 
                              ovx_choice_rows(), 
                              region_choice_rows(), 
                              exclude_samples_rows()))
    meta_trap[final_rows]
  })
  
  #Provide comparison options and create a filter based on whether there is more than one factor left in the filtered metadata
  comparison_options <- c("age", "ovx", "region", "sex")
  comparison_options_filter <- reactive({
    lapply(meta()[, ..comparison_options], function(x) length(unique(x))) > 1
  })
  #Provide covariate options and create a filter based on whether there is more than one factor left in the filtered metadata
  covariates_options <- c("th_enrichment", "age", "sex", "ovx", "region")
  covariates_options_filter <- reactive({
    lapply(meta()[, ..covariates_options], function(x) length(unique(x))) > 1
  })

  comparison_choices_listen <- reactive({
    list(input$trap_choice,
         input$remove_outliers,
         input$age_choice,
         input$ovx_choice, 
         input$region_choice, 
         input$exclude_samples_choice)
  })

  #Update the comparison choices available based on the TRAP experiment selected
  observeEvent(comparison_choices_listen(), {
    updateSelectizeInput(session, 
                         "comparison_choice", 
                         choices = {
                           options <- comparison_options[comparison_options_filter()]
                           }, 
                         selected = lapply(reactiveValuesToList(input), unclass)$comparison_choice)
  })
  
  nominator_choices_listen <- reactive({
    list(input$trap_choice,
         input$remove_outliers,
         input$age_choice,
         input$ovx_choice, 
         input$region_choice, 
         input$exclude_samples_choice, 
         input$comparison_choice)
  })
  
  #Update the nominator choices available based on the comparison selected
  observeEvent(nominator_choices_listen(), {
    updateSelectizeInput(session,
                         "nominator_choice",
                         choices = unique(meta()[[input$comparison_choice]]), 
    selected = lapply(reactiveValuesToList(input), unclass)$nominator_choice)
    })
  
  #Update the demoninator choices available based on the nominator selected
  observeEvent(input$nominator_choice, {
    updateSelectizeInput(session, 
                         "denominator_choice", 
                         choices = {
                           options <- unique(meta()[[input$comparison_choice]])
                           options[!options == input$nominator_choice]
                         }, 
                         selected = lapply(reactiveValuesToList(input), unclass)$denominator_choice)
  })
  

  
  #This listening function allow the covariates options to update when any of the below listed inputs changes. 
  # I could add the other options (age, region), but it might be better to write a function that checks for 
  # errors at the end instead and stops the button actioning.
  covariates_choices_listen <- reactive({
    list(input$trap_choice,
         input$remove_outliers,
         input$age_choice,
         input$ovx_choice, 
         input$region_choice, 
         input$exclude_samples_choice,
         input$comparison_choice)
  })
  
  observeEvent(covariates_choices_listen(), {
    updateSelectizeInput(session, 
                         "covariates_choice", 
                         choices = {
                           options <- covariates_options[covariates_options_filter()]
                           options[options != input$comparison_choice] 
                         },
                         selected = "th_enrichment")
  })
  
  # Filter data based on selections
  output$meta_select <- DT::renderDataTable({
    meta()
  })

  output$ui_makedds <- renderUI({
    actionButton("button_makedds", "Generate the dataset", class = "btn btn-success")
  })
  
  generate_dds <- reactive({
    #baseMean for filtering before independent filtering
    #This influences WGCNA results, because of the resulting matrix size
    basemean <- 10
    
    files <- file.path("/home/ubuntu/isitdifferentiallyexpressed/input_data/kallisto", meta()$code, "abundance.h5")
    names(files) <- meta()$name
    txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = T)
    dds <- DESeqDataSetFromTximport(txi, colData = meta(), design = as.formula(paste0("~", paste(input$covariates_choice, input$comparison_choice, sep = " + "), collapse = " + ")))
    dds <- estimateSizeFactors(dds)
    keep_feature <- rowMeans(counts(dds, normalized = TRUE)) >= basemean
    dds_removed <- dds[!keep_feature, ]
    dds <- dds[keep_feature, ]
    dds <- DESeq(dds,
                 minReplicatesForReplace = Inf,
                 parallel = TRUE)

    return(dds)
  })
  
  generate_counts <- reactive({
    
  })
  
  generate_res <- reactive({
    # Set FDR
    fdr = 0.05
    # Set lfc at 5%
    lfc = log2(1.05)
    #Filtering with independent hypothesis weighting
    filterfunction <- "ihw"
    #Ensuring independent filtering is performed
    independentfiltering = TRUE
    #Using apeglm shrinkage (note this requires use of the coef argument as opposed to contrast for interaction results (see ?results))
    shrinkage = "apeglm"
    
    res <- results(values$dds,
                   contrast = c(input$comparison_choice, input$nominator_choice, input$denominator_choice),
                   alpha = fdr,
                   lfcThreshold = lfc,
                   filterFun = get(filterfunction),
                   independentFiltering = independentfiltering,
                   parallel = TRUE)
    return(res)
  })
  
  generate_res_summary <- reactive({
    # Set FDR
    res_summary <- make_res_summary(values$res, 0.05)
    return(res_summary)
  })
  
  observeEvent(input$button_makedds,
               {
                 values$dds <- generate_dds()
                 values$res <- generate_res()
                 values$res_summary <- generate_res_summary()
               })
  
  output$res_summary <- DT::renderDataTable({
    if(!is.null(values$res_summary)){
      values$res_summary
    }
  })
  
  # output$results_plot <- renderPlotly({
  #   if(!is.null(values$res_summary)){
  #     g <- make_tidy(values$res, "gene") %>%
  #     ggplot(aes(y = -log10(padj), 
  #                x = log2FoldChange, 
  #                label = gene)) +
  #     geom_point(alpha = 0.25, 
  #                size = 2) +
  #       theme_cowplot() +
  #       scale_color_lancet()
  #     # return(g)
  #     ggplotly(g)
  #   }
  # })
  
  output$results_plot <- renderPlotly({
    if(!is.null(values$res_summary)){
      plot_ly(data = make_tidy(values$res, "gene") %>%
                mutate(padj = -log10(padj))) %>% 
        add_trace(x = ~padj, 
                  y = ~log2FoldChange,
                  type = "scattergl", 
                  mode = "markers", 
                  alpha = 0.3)
      }
  })

  output$debugdiy <- renderPrint({
    if(!is.null(values$res)){
      # print(values$dds_obj)
      # print(design(values$dds_obj))
      # print(resultsNames(values$dds_obj))
      # print(input$comparison_choice)
      # print(input$nominator_choice)
      # print(input$denominator_choice)
      print(summary(values$res))
    }
  })
  
  
  # output$dds_design <- renderUI({
  #   if(is.null(values$exp_design))
  #     return(NULL)
  #   covars <- colnames(values$exp_design)
  #   selectInput("dds_design", 
  #               label = "Select covariates to model: ", 
  #               choices = c(NULL, covars), 
  #               selected = NULL, 
  #               multiple = TRUE)
  # })
  
}

shinyApp(ui, server)