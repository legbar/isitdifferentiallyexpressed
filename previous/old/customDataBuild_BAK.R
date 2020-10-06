library(tximport)
library(DESeq2)
library(rhdf5)

meta_trap <- readRDS("data/meta_trap.rds")
tx2gene_trap <- readRDS("data/tx2gene_trap.rds")

customDataBuildUI <- function(id){
  ns <- NS(id)
  tagList(
    selectizeInput(
      inputId = ns("dataset"), 
      label = "Choose a Dataset", 
      choices = c("TRAP", "iPSC"), 
      selected = NULL),
    conditionalPanel(
      ns = ns,
      condition = 'input.dataset == "TRAP"',
      customDataBuildTRAPUI("customDataBuildTRAP"),
    ),
    # conditionalPanel(
    #   ns = ns,
    #   condition = 'input.dataset == "TRAP"',
    #   selectizeInput(
    #     inputId = ns("trapCohortSelect"),
    #     label = "Which TRAP Cohort?",
    #     choices = unique(meta_trap$cohort)
    #     )
    #   ),
    conditionalPanel(
      ns = ns,
      condition = 'input.dataset == "TRAP"',
    checkboxGroupInput(
      inputId = ns("ageSelect"),
      label = "Which age groups?",
      choices = unique(meta_trap$Age),
      select = unique(meta_trap$Age)
    )
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

    selectizeInput(
      inputId = ns("excludeSelect"),
      label = "Exclude any samples?",
      choices = NULL,
      multiple = TRUE),
    selectizeInput(
      inputId = ns("comparisonSelect"),
      label = "Comparison Factor",
      choices = NULL),
    selectizeInput(
      inputId = ns("covariateSelect"),
      label = "Covariates",
      choices = NULL, 
      multiple = TRUE),
    actionButton(
      inputId = ns("buttonMakeDDS"), 
      label = "Generate the dataset"),
    hr(),
    verbatimTextOutput(
      ns("debugdiy")
      )
  )
}

customDataBuildServer <- function(id){
  moduleServer(
    id,
    function(input, output, session){
      
      values <- reactiveValues()
      
      values$trapCohortRows <- customDataBuildTRAPServer("customDatasetBuildTRAP")
      
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
      
      values$metaRowFilter <- reactive({
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
        meta_trap[rows,]
      })
      
      observeEvent(values$metaRowFilter(), {
        updateSelectizeInput(session,
                             "excludeSelect",
                             choices = values$metaRowFilter()$name,
                             selected = lapply(reactiveValuesToList(input), unclass)$excludeSelect)
      })
      
      #Rows excluded specifically
      values$includeRows <- reactive({
        which(!meta_trap$name %in% input$excludeSelect)
      })


      #Build the filtered table by intersecting a list of all the row numbers
      values$metaFiltered <- reactive({
        rows <- Reduce(intersect,
                       list(values$trapCohortRows(),
                            values$outlierRows(),
                            values$ipRows(),
                            values$ageRows(),
                            values$ovxRows(),
                            values$includeRows()
                       )
        )
        meta_trap[rows,]
      })
      
      #Provide comparison options and create a filter based on whether there is more than one factor left in the filtered metadata
      values$comparisonOptions <- reactive({
        options <- c("Age", "Genotype", "Region", "Sex")
        filter <- lapply(values$metaFiltered()[, options], function(x) length(unique(x))) > 1
        options[filter]
      })
      
      #Update the comparison choices available based on the TRAP experiment selected
      observeEvent(values$metaFiltered(), {
        updateSelectizeInput(session,
                             "comparisonSelect",
                             choices = values$comparisonOptions(),
                             selected = NULL)
      })
     
      #Provide covariate options and create a filter based on whether there is more than one factor left in the filtered metadata
      values$covariateOptions <- reactive({
        options <- c("Age", "Genotype", "Region", "Sex", "th_enrichment")
        filter <- lapply(values$metaFiltered()[, options], function(x) length(unique(x))) > 1
        options[filter]
      })

      observeEvent(input$denominatorSelect, {
        updateSelectizeInput(session,
                             "covariateSelect",
                             choices = {
                               options <- values$covariateOptions()
                               options[options != input$comparisonSelect]
                             },
                             selected = "th_enrichment")
      })

      generateDDS <- reactive({
        #baseMean for filtering before independent filtering
        #This influences WGCNA results, because of the resulting matrix size
        basemean <- 10          
        files <- file.path("/zfs/analysis/kallisto_store", values$metaFiltered()$code, "abundance.h5")
        names(files) <- values$metaFiltered()$name
        txi <- tximport(files, type = "kallisto", tx2gene = tx2gene_trap, ignoreTxVersion = T)

        dds <- DESeqDataSetFromTximport(txi, colData = values$metaFiltered(), design = as.formula(paste0("~", paste(input$covariateSelect, input$comparisonSelect, sep = " + "), collapse = " + ")))
        # dds[[input$comparisonSelect]] <- relevel(dds[[input$comparisonSelect]], input$denominatorSelect)
        dds <- estimateSizeFactors(dds)
        keep_feature <- rowMeans(counts(dds, normalized = TRUE)) >= basemean
        dds_removed <- dds[!keep_feature, ]
        dds <- dds[keep_feature, ]
        dds <- DESeq(dds,
                     minReplicatesForReplace = Inf,
                     parallel = TRUE)

        print(resultsNames(dds))

        return(dds)
      })
      
      observeEvent(input$buttonMakeDDS,
                   {
                     values$dds <- generateDDS()
                   })
      
      output$debugdiy <- renderPrint({
        if(!is.null(values$dds)){
          print(values$dds)
          print(resultsNames(values$dds))
          print(input$comparisonSelect)

        }
      })
      
      # observeEvent(input$buttonMakeDDS, {
      #   updateTabsetPanel(session, 
      #                     "customDatasetHiddenTabs", 
      #                     selected = "View Results")
      # })
      
      return(values)
      }
  )
}

