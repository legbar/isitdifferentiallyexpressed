library(tximport)
library(DESeq2)
library(rhdf5)

customDataBuildGenericUI <- function(id) {
  ns <- NS(id)
  tagList(
    selectizeInput(
      inputId = ns("experimentSelect"), 
      label = "Select an Experiment Type", 
      choices = c("TRAP", "iPSC"), 
      selected = "iPSC"
    ),
    conditionalPanel(
      ns = ns,
      condition = 'input.experimentSelect == "TRAP"',
      customDataBuildSpecificUI("customDataBuild")
    ),
    selectizeInput(
      inputId = ns("excludeSelect"),
      label = "Exclude any samples?",
      choices = NULL,
      multiple = TRUE
    ),
    selectizeInput(
      inputId = ns("comparisonSelect"),
      label = "Comparison Category",
      choices = NULL
    ),
    selectizeInput(
      inputId = ns("covariateSelect"),
      label = "Covariates",
      choices = NULL,
      multiple = TRUE
    ),
    actionButton(inputId = ns("buttonMakeDDS"),
                 label = "Generate the dataset"),
    verbatimTextOutput(ns("debug"))
  )
}


customDataBuildGenericServer <- function(id) {
  moduleServer(id,
               function(input, output, session) {
                 
                 values <- customDataBuildSpecificServer("customDataBuild")
                 
                 observeEvent(values$metaSpecificRowFilter(), {
                   updateSelectizeInput(
                     session,
                     "excludeSelect",
                     choices = values$metaSpecificFiltered()$name,
                     selected = NULL
                   )
                 })
                 
                 #Rows excluded specifically
                 values$includeRows <- reactive({
                   which(!values$meta()$name %in% input$excludeSelect)
                 })
                 
                 #Build the filtered table by intersecting a list of all the row numbers
                 values$metaFiltered <- reactive({
                   rows <- Reduce(intersect,
                                  list(values$metaSpecificRowFilter(),
                                       values$includeRows()))
                   values$meta()[rows, ]
                 })
                 
                 #Update the comparison choices available based on the TRAP experiment selected
                 observeEvent(values$metaFiltered(), {
                   updateSelectizeInput(
                     session,
                     "comparisonSelect",
                     choices = values$comparisonOptions(),
                     selected = NULL
                   )
                 })
                 
                 observeEvent(input$comparisonSelect, {
                   updateSelectizeInput(session,
                                        "covariateSelect",
                                        choices = {
                                          options <- values$covariateOptions()
                                          options[options != input$comparisonSelect]
                                        },
                                        selected = "th_enrichment")
                 })
                 
                 values$comparisonSelect <- reactive(input$comparisonSelect)
                 
                 generateDDS <- reactive({
                   #baseMean for filtering before independent filtering
                   #This influences WGCNA results, because of the resulting matrix size
                   basemean <- 10
                   files <-
                     file.path("/zfs/analysis/kallisto_store",
                               values$metaFiltered()$code,
                               "abundance.h5")
                   names(files) <- values$metaFiltered()$name
                   txi <-
                     tximport(
                       files,
                       type = "kallisto",
                       tx2gene = values$tx2gene(),
                       ignoreTxVersion = T
                     )
                   
                   dds <-
                     DESeqDataSetFromTximport(txi,
                                              colData = values$metaFiltered(),
                                              design = as.formula(paste0(
                                                "~",
                                                paste(input$covariateSelect, input$comparisonSelect, sep = " + "),
                                                collapse = " + "
                                              )))
                   # dds[[input$comparisonSelect]] <- relevel(dds[[input$comparisonSelect]], input$denominatorSelect)
                   dds <- estimateSizeFactors(dds)
                   keep_feature <-
                     rowMeans(counts(dds, normalized = TRUE)) >= basemean
                   dds_removed <- dds[!keep_feature,]
                   dds <- dds[keep_feature,]
                   dds <- DESeq(dds,
                                minReplicatesForReplace = Inf,
                                parallel = TRUE)
                   print("DDS Object build complete.")
                   return(dds)
                 })
                 
                 observeEvent(input$buttonMakeDDS,
                              {
                                values$dds <- generateDDS()
                              })
                 
                 output$debug <- renderPrint({
                   
                 })
                 
                 return(values)
               })
}
