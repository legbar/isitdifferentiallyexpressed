library(shiny)
library(tidyverse)
library(shinydashboard)
library(data.table)


meta_trap <- readRDS("data/meta_trap.rds")
setDT(meta_trap)
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
                                    choices = unique(meta_trap$name), 
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
                     verbatimTextOutput("debugdiy")
                     ),
              column(
                width = 9, 
                box(
                  title = "TRAP Metadata",
                  DT::dataTableOutput("meta_select"), 
                  width = NULL, 
                  solidHeader = TRUE
                  )
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
  
  values <- reactiveValues()
  
  trap_choice_rows <- reactive({
    meta_trap[cohort %in% input$trap_choice, which = TRUE]
  })
  outlier_rows <- reactive({
    if (input$remove_outliers == TRUE) {
    meta_trap[outlier == FALSE, which = TRUE]
    } else {
    meta_trap[outlier == TRUE | outlier == FALSE, which = TRUE]
    }
  })
  ip_rows <- reactive({
    if (input$ip_only == TRUE) {
      meta_trap[ip == "ip", which = TRUE]
    } else {
      meta_trap[ip %in% meta_trap$ip, which = TRUE]
    }
  })
  age_choice_rows <- reactive({
    meta_trap[age %in% input$age_choice, which = TRUE]
  })
  ovx_choice_rows <- reactive({
    meta_trap[ovx %in% input$ovx_choice, which = TRUE]
  })
  region_choice_rows <- reactive({
    meta_trap[region %in% input$region_choice, which = TRUE]
  })
  exclude_samples_rows <- reactive({
    meta_trap[!name %in% input$exclude_samples_choice, which = TRUE ]
  })
  
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
  
  comparison_options <- c("age", "ovx", "region", "sex")
  comparison_options_filter <- reactive({
    lapply(meta()[, ..comparison_options], function(x) length(unique(x))) > 1
  })
  
  covariates_options <- c("th_enrichment", "age", "sex", "ovx", "region")
  covariates_options_filter <- reactive({
    lapply(meta()[, ..covariates_options], function(x) length(unique(x))) > 1
  })
  
  observeEvent(input$trap_choice, {
    updateSelectizeInput(session, 
                         "comparison_choice", 
                         choices = comparison_options[comparison_options_filter()])
  })
    
  observeEvent(input$age_choice, {
    updateSelectizeInput(session, 
                         "comparison_choice", 
                         choices = comparison_options[comparison_options_filter()])
  })
  
  observeEvent(input$ovx_choice, {
    updateSelectizeInput(session, 
                         "comparison_choice", 
                         choices = comparison_options[comparison_options_filter()])
  })
  
  observeEvent(input$region_choice, {
    updateSelectizeInput(session, 
                         "comparison_choice", 
                         choices = comparison_options[comparison_options_filter()])
  })
  
  observeEvent(input$comparison_choice, {
    updateSelectizeInput(session, 
                         "nominator_choice", 
                         choices = unique(meta()[[input$comparison_choice]]))
    })
  
  observeEvent(input$nominator_choice, {
    updateSelectizeInput(session, 
                         "denominator_choice", 
                         choices = {
                           options <- unique(meta()[[input$comparison_choice]])
                           options[!options == input$nominator_choice]
                         })
  })
  
  #This listening function allow the covariates options to update when any of the below listed inputs changes. 
  # I could add the other options (age, region), but it might be better to write a function that checks for 
  # errors at the end instead and stops the button actioning.
  covariates_listen <- reactive({
    list(input$trap_choice, 
         input$comparison_choice, 
         input$ovx_choice)
  })
  
  observeEvent(covariates_listen(), {
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
  
  diyDDS <- reactive({

    files <- file.path("/zfs/analysis/thesis-pk/input_data/kallisto_PK1_KW2/", meta()$code, "abundance.h5")
    names(files) <- meta()$name
    txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = T)
    dds <- DESeqDataSetFromTximport(txi, colData = meta(), design = as.formula(paste0("~", paste(input$covariates_choice, input$comparison_choice, sep = " + "), collapse = " + ")))
    # dds <- DESeqDataSetFromMatrix(countData = values$countmatrix,
    #                               colData = values$expdesign,
    #                               design=as.formula(paste0("~",paste(input$dds_design, collapse=" + "))))

    dds <- estimateSizeFactors(dds)
    return(dds)
  })

  observeEvent(input$button_makedds,
               {
                 values$dds_obj <- diyDDS()
               })

  output$debugdiy <- renderPrint({
    if(!is.null(values$dds_obj)){
      print(values$dds_obj)
      print(design(values$dds_obj))
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