server <- function(input, output, session) {
  datasetInput <- reactive({
    switch(input$dataset, 
           "TRAP" = trap, 
           "iPSC")
  })
  
  #Make a list of reactive values - not currently in use
  values <- reactiveValues()