meta_trap <- readRDS("data/meta_trap.rds") %>%
  filter(outlier == FALSE)
tx2gene_trap <- readRDS("data/tx2gene_trap.rds")

# TRAP DDS Build Questions
ddsBuild_TRAP_UI <- function(id) {
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

    checkboxGroupInput(
      inputId = ns("ageSelect"),
      label = "Which age groups?",
      choices = unique(meta_trap$Age),
      select = unique(meta_trap$Age)
    )
  )
}


#Server functions
ddsBuild_TRAP_Server <- function(id) {
  moduleServer(id,
               function(input, output, session) {
                 
                 values <- reactiveValues()
                 
                 #Cohort rows
                 values$trapCohortRows <- reactive({
                   which(meta_trap$cohort %in% input$trapCohortSelect)
                 })
                 
                #Age rows
                 values$ageRows <- reactive({
                   which(meta_trap$Age %in% input$ageSelect)
                 })
                 
                 #Genotype rows
                 values$ovxRows <- reactive({
                   which(meta_trap$Genotype %in% input$ovxSelect)
                 })
                 
                 #Region rows
                 values$regionRows <- reactive({
                   which(meta_trap$Region %in% input$regionSelect)
                 })
                 
                 #Filter based on rows
                 values$rowFilter <- reactive({
                 rows <- Reduce(
                   intersect,
                   list(
                     values$trapCohortRows(),
                     values$ageRows(),
                     values$ovxRows(),
                     values$regionRows()
                   )
                 )
                 })
                 
                 # filter meta
                  values$metaFiltered <- reactive({
                    meta_trap[values$rowFilter(),]
                  })
                  
                  #Provide comparison options and create a filter based on whether there is more than one factor left in the filtered metadata
                  values$comparisonOptions <- reactive({
                    options <- c("Age" = "age", "Genotype" = "ovx", "Region" = "region", "Sex" = "sex")
                    filter <-
                      lapply(values$metaFiltered()[, options], function(x)
                        length(unique(x))) > 1
                    options[filter]
                  })
                  
                  #Provide covariate options and create a filter based on whether there is more than one factor left in the filtered metadata
                  values$covariateOptions <- reactive({
                    options <- c("Age" = "age", "Genotype" = "ovx", "Region" = "region", "Sex" = "sex", "TH Enrichment" = "th_enrichment")
                    filter <-
                      lapply(values$metaFiltered()[, options], function(x)
                        length(unique(x))) > 1
                    options[filter]
                  })
                  
                  values$displayColumns <- reactive(c("name",
                                             "cohort", 
                                             "Age",
                                             "Region",
                                             "Genotype"))
                  
                  return(values)
  
               }
)
}




  
  
  
