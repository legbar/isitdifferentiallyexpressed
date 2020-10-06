# DDS Build Skeleton
ddsBuild_skeleton_UI <-
  function(id) {
    ns <- NS(id)
    
    tagList(
      selectizeInput(
        inputId = ns("experimentSelect"),
        label = "Select an Experiment Type",
        choices = c("TRAP", "iPSC"),
        selected = "TRAP"
      ),
      
      uiOutput(ns("experimentSpecificUI")),
      
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
                   label = "Generate the dataset")
      
    )
    
  }

#Server functions
ddsBuild_skeleton_Server <- function(id) {
  moduleServer(id,
               function(input, output, session) {
                 
                 output$experimentSpecificUI <- renderUI({
                   ns <- session$ns
                   req(input$experimentSelect)
                   if (input$experimentSelect == "TRAP"){
                     ddsBuild_TRAP_UI(ns("ddsBuild_TRAP"))
                   } else {
                     ddsBuild_iPSC_UI(ns("ddsBuild_iPSC"))
                   }
                 })
                 
                 values <- ddsBuild_TRAP_Server("ddsBuild_TRAP")

                 observeEvent(values$metaFiltered(), {
                   updateSelectizeInput(
                     session,
                     "excludeSelect",
                     choices = values$metaFiltered()$name,
                     selected = NULL
                   )
                 })

                 #Rows excluded specifically
                 values$includeRows <- reactive({
                   which(!values$metaFiltered()$name %in% input$excludeSelect)
                 })

                 #Build the filtered table by intersecting a list of all the row numbers
                 values$metaFinal <- reactive({
                   values$metaFiltered()[values$includeRows(), ]
                 })

                 #Update the comparison choices available based on the TRAP experiment selected
                 observe({
                   updateSelectizeInput(
                     session,
                     "comparisonSelect",
                     choices = values$comparisonOptions(),
                     selected = NULL
                   )
                 })

                 observe({
                   updateSelectizeInput(session,
                                        "covariateSelect",
                                        choices = {
                                          options <- values$covariateOptions()
                                          options[options != input$comparisonSelect]
                                        },
                                        selected = "th_enrichment")
                 })
                 
                 values$covariateSelect <- reactive(input$covariateSelect)
                 values$comparisonSelect <- reactive(input$comparisonSelect)
                 
                 values$comparisonType <- reactive({
                   if (is.factor(values$metaFinal()[[values$comparisonSelect()]])){
                     "FACTOR"
                   } else {
                     "NUMERIC"
                   }
                 })
                 
                 values$dds <- eventReactive(input$buttonMakeDDS,
                              {
                                withProgress(message = 'Building Dataset', value = 0, {
                                generateDDS(files = file.path("/zfs/analysis/kallisto_store", 
                                                                            values$metaFinal()$code, 
                                                                            "abundance.h5"), 
                                                          meta = values$metaFinal(), 
                                                          TX2GENE = tx2gene_trap, 
                                                          quantType = "kallisto", 
                                                          covariates = values$covariateSelect(), 
                                                          comparison = values$comparisonSelect(), 
                                                          basemean = 10)
                              })
                              })
                 
                 values$counts <- eventReactive(input$buttonMakeDDS,
                                                {
                                                  counts(values$dds(), normalized = T) %>%
                                                    as_tibble(rownames = "gene") %>%
                                                    pivot_longer(-gene, names_to = "name", values_to = "count") %>%
                                                    inner_join(meta_trap, by = "name")
                                                })
                 
                 values$vsd <- eventReactive(input$buttonMakeDDS,
                                                {
                                                  vst(values$dds(), blind = FALSE) 
                                                })
                 
                 values$vsd_long_meta <- reactive({
                   values$vsd() %>%
                     assay(.) %>%
                     as_tibble(rownames = "gene") %>%
                     pivot_longer(-gene, names_to = "name", values_to = "count") %>%
                     inner_join(meta_trap, by = "name")
                 })
                 
                 values$pca <- reactive({
                   vsd_var <- assay(values$vsd())[rowVars(assay(values$vsd())) > 0,]
                   pca <- prcomp(t(vsd_var)) #Not scaling as passing log transformed data: https://www.biostars.org/p/280615/#283126
                  })
                 
                 values$pca_meta <- reactive({
                   pca_meta <- values$pca()$x %>%
                     as_tibble(rownames = "name") %>%
                     inner_join(values$metaFinal(), by = "name") 
                 })
                 
                 values$pca_percentVar <- reactive({
                   values$pca()$sdev^2/sum(values$pca()$sdev^2) #Calculate PC variance
                 })
                 
                   

                 return(values)
                 
               })
}