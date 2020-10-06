customDataViewUI <- function(id, label = "Custom Dataset View"){
  ns <- NS(id)
  
  fluidRow(
    column(width = 3, 
           box(
             title = "Introduction", 
             height = NULL, 
             solidHeader = TRUE, 
             status = "info", 
             width = NULL, 
             p("Use the options below to select a dataset and compare groups of interest.")
           )
           ), 
    column(
      width = 9,
      box(
        title = "TRAP Metadata",
        DT::dataTableOutput(ns("metaFiltered")),
        width = NULL,
        solidHeader = TRUE
      )
    ) 
  )
}

customDataViewServer <- function(id){
  moduleServer(
    id,
    function(input, output, session){
      
      # Filter data based on selections
      output$metaFiltered <- DT::renderDataTable({
        values$metaFiltered() %>%
          select("Name" = name, 
                 "Cohort" = cohort,
                 Region,
                 Age,
                 Genotype)
      })
      
    }
  )
}
