
customDatasetUI <- function(id, label = "Custom Dataset"){
  # `NS(id)` returns a namespace function, which was save as `ns` and will
  # invoke later.
  ns <- NS(id)

  
          fluidRow(
            column(width = 3, 
                   box(title = "Introduction", 
                       height = NULL, 
                       solidHeader = TRUE, 
                       status = "info", 
                       width = NULL, 
                       p("Use the options below to select a dataset and compare groups of interest.")
                   ), 
                   selectInput(inputId = ns("dataset"), 
                               label = "Choose a Dataset", 
                               choices = c("TRAP"), 
                               selected = NULL,
                               selectize = TRUE), 
                   #TRAP choices
                   conditionalPanel(
                     condition = 'input.dataset == "TRAP"', 
                     selectizeInput(inputId = ns("trapCohort"), 
                                    label = "Which TRAP Cohort?", 
                                    choices = unique(meta_trap$cohort))),
                   conditionalPanel(
                     condition = 'input.dataset == "TRAP"', 
                     checkboxGroupInput(inputId = ns("age_choice"), 
                                        label = "Which age groups?", 
                                        choices = unique(meta_trap$age), 
                                        select = unique(meta_trap$age))),
                   conditionalPanel(
                     condition = 'input.trapCohort != "PK1"', 
                     checkboxGroupInput(inputId = ns("ovx_choice"), 
                                        label = "Which genotypes?", 
                                        choices = unique(meta_trap$ovx), 
                                        select = unique(meta_trap$ovx))),
                   conditionalPanel(
                     condition = 'input.dataset == "TRAP" & input.trapCohort == "PK1"', 
                     checkboxGroupInput(inputId = ns("region_choice"), 
                                        label = "Which regions?", 
                                        choices = unique(meta_trap$region), 
                                        select = unique(meta_trap$region))),
                   conditionalPanel(
                     condition = 'input.dataset == "TRAP"', 
                     checkboxInput(inputId = ns("ip_only"),
                                   label = "TRAP Samples Only", 
                                   value = TRUE)),
                   conditionalPanel(
                     condition = 'input.trapCohort == "PK1"', 
                     checkboxInput(inputId = ns("remove_outliers"), 
                                   label = "Remove outliers", 
                                   value = TRUE)),
                   selectizeInput(inputId = ns("exclude_samples_choice"), 
                                  label = "Exclude any samples?", 
                                  choices = NULL, 
                                  multiple = TRUE), 
                   selectizeInput(inputId = ns("comparison_choice"), 
                                  label = "Comparison Factor", 
                                  choices = NULL),
                   fluidRow(
                     column(
                       width = 4, 
                       selectizeInput(inputId = ns("nominator_choice"), 
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
                       selectizeInput(inputId = ns("denominator_choice"), 
                                      label = "Denominator", 
                                      choices = NULL)
                     )
                   ),
                   selectizeInput(inputId = ns("covariates_choice"), 
                                  label = "Covariates to include:", 
                                  multiple = TRUE,
                                  choices = NULL),
                   verbatimTextOutput(ns("final_rows_value")),
                   uiOutput(ns("ui_makedds")), 
                   hr(),
                   verbatimTextOutput(ns("debugdiy"))
                   # box(
                   #   title = "Result Summary", 
                   #   DT::dataTableOutput("res_summary"), 
                   #   width = NULL, 
                   #   solidHeader = TRUE
                   # )
            ),
            column(
              width = 9, 
              box(
                title = "TRAP Metadata",
                DT::dataTableOutput(ns("meta_select")), 
                width = NULL, 
                solidHeader = TRUE
              )
              # box(title = "Results Plot", 
              #     plotlyOutput("results_plot"))
              # )
            )
          )

}

customDatasetServer <- function(id){
  moduleServer(
    id,
    function(input, output, session){
      
      #Make a list of reactive values - not currently in use
      values <- reactiveValues()
      
      #Rows selected by TRAP choice
      values$trapCohortRows <-  reactive({
        meta_trap[cohort %in% input$trapCohort, which = TRUE]
      })
      
      #Rows selected by excluding outliers
      values$outlier_rows <- reactive({
        if (input$remove_outliers == TRUE) {
          meta_trap[outlier == FALSE, which = TRUE]
        } else {
          meta_trap[outlier == TRUE | outlier == FALSE, which = TRUE]
        }
      })
      #Rows selected by choosing only IP samples
      values$ip_rows <- reactive({
        if (input$ip_only == TRUE) {
          meta_trap[ip == "ip", which = TRUE]
        } else {
          meta_trap[ip %in% meta_trap$ip, which = TRUE]
        }
      })
      #Rows from age_choice
      values$age_rows <- reactive({
        meta_trap[age %in% input$age_choice, which = TRUE]
      })
      #Rows from OVX choice
      values$ovx_rows <- reactive({
        meta_trap[ovx %in% input$ovx_choice, which = TRUE]
      })
      #Rows from region choice
      values$region_rows <- reactive({
        meta_trap[region %in% input$region_choice, which = TRUE]
      })
      
      meta0 <- reactive({
        rows <- Reduce(intersect, 
                       list(values$trapCohortRows(),
                            values$outlier_rows(),
                            values$ip_rows(), 
                            values$age_rows(), 
                            values$ovx_rows(), 
                            values$region_rows()))
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
                             list(values$trapCohortRows(),
                                  values$outlier_rows(),
                                  values$ip_rows(), 
                                  values$age_rows(), 
                                  values$ovx_rows(), 
                                  values$region_rows(), 
                                  exclude_samples_rows()))
        meta_trap[final_rows]
      })
      
    }
  )
}

