meta_trap <- readRDS("data/meta_trap.rds")
tx2gene_trap <- readRDS("data/tx2gene_trap.rds")

customDataBuildSpecificUI <- function(id){
  ns <- NS(id)
  
  tagList(
    
        selectizeInput(
          inputId = ns("trapCohortSelect"),
          label = "Which TRAP Cohort?",
          choices = unique(meta_trap$cohort)
        ),
    
    conditionalPanel(
      ns = ns,
      condition = 'input.trapCohortSelect != "PK1"',
      checkboxGroupInput(
        inputId = ns("ovxSelect"),
        label = "Which genotypes?",
        choices = unique(meta_trap$Genotype),
        select = unique(meta_trap$Genotype)
      )
    ),
    
    conditionalPanel(
      ns = ns,
      condition = 'input.trapCohortSelect == "PK1"',
      checkboxGroupInput(
        inputId = ns("regionSelect"),
        label = "Which regions?",
        choices = unique(meta_trap$Region),
        select = unique(meta_trap$Region)
      )
    ),
    
    conditionalPanel(
      ns = ns,
      condition = 'input.dataset == "TRAP"',
      checkboxInput(
        inputId = ns("ipSelect"),
        label = "TRAP Samples Only",
        value = TRUE
      )
    ),
    
    conditionalPanel(
      ns = ns,
      condition = 'input.trapCohortSelect == "PK1"',
      checkboxInput(
        inputId = ns("outlierSelect"),
        label = "Include outliers",
        value = FALSE
      )
    ),
    
    checkboxGroupInput(
      inputId = ns("ageSelect"),
      label = "Which age groups?",
      choices = unique(meta_trap$Age),
      select = unique(meta_trap$Age)
    ),
    
    hr()
    )
}
  
customDataBuildSpecificServer <- function(id){
    moduleServer(
      id,
      function(input, output, session){
        
        values <- reactiveValues()
        
        values$trapCohortRows <- reactive({
          which(meta_trap$cohort %in% input$trapCohortSelect)
        })
        
        #Rows selected by excluding outliers
        values$outlierRows <- reactive({
          if (input$outlierSelect == TRUE) {
            seq(1, nrow(meta_trap), by = 1)
          } else {
            which(meta_trap$outlier == FALSE)
          }
        })
        
        #Rows selected by choosing only IP samples
        values$ipRows <- reactive({
          if (input$ipSelect == TRUE) {
            which(meta_trap$ip == "ip")
          } else {
            seq(1, nrow(meta_trap), by = 1)
          }
        })
        
        #Rows from age_choice
        values$ageRows <- reactive({
          which(meta_trap$Age %in% input$ageSelect)
        })
        
        #Rows from OVX choice
        values$ovxRows <- reactive({
          which(meta_trap$Genotype %in% input$ovxSelect)
        })
        
        #Rows from region choice
        values$regionRows <- reactive({
          which(meta_trap$Region %in% input$regionSelect)
        })
        
        values$metaSpecificRowFilter <- reactive({
          rows <- Reduce(intersect,
                         list(
                           values$trapCohortRows(),
                           values$outlierRows(),
                           values$ipRows(),
                           values$ageRows(),
                           values$ovxRows(),
                           values$regionRows()
                         )
          )
        })
        
        values$metaSpecificFiltered <- reactive({
          meta_trap[values$metaSpecificRowFilter(),]
          })
        
        #Provide comparison options and create a filter based on whether there is more than one factor left in the filtered metadata
        values$comparisonOptions <- reactive({
          options <- c("Age", "Genotype", "Region", "Sex")
          filter <- lapply(values$metaFiltered()[, options], function(x) length(unique(x))) > 1
          options[filter]
        })
        
        #Provide covariate options and create a filter based on whether there is more than one factor left in the filtered metadata
        values$covariateOptions <- reactive({
          options <- c("Age", "Genotype", "Region", "Sex", "th_enrichment")
          filter <- lapply(values$metaFiltered()[, options], function(x) length(unique(x))) > 1
          options[filter]
        })
        
        return(values)
        
      }
    )
  }
  