customDataViewUI <- function(id){
  ns <- NS(id)
  fluidRow(
    column(
      width = 4, 
      box(
        title = "Configure Results",
        selectizeInput(
          inputId = ns("contrastSelect"),
          label = "Select specific comparison",
          choices = NULL,
          multiple = FALSE),
        width = TRUE,
        solidHeader = TRUE
      ), 
      verbatimTextOutput(
        ns("debug")
      )
    )
  )
  
}

customDataViewServer <- function(id, buildValues){
  moduleServer(
    id,
    function(input, output, session){
      
      values <- buildValues
      
      values$comb <- reactive({
        combinations <- c(combn(unique(values$metaFiltered()[[values$comparisonSelect()]]), 2, simplify = F),
                                    lapply(combn(unique(values$metaFiltered()[[values$comparisonSelect()]]), 2, simplify = F), function(x){rev(x)}))
        
        contrasts <- lapply(combinations, function(x){paste(x, collapse = " versus ")})
      })
      
      observeEvent(values$dds, {
        updateSelectizeInput(session,
                             "contrastSelect",
                             choices = values$comb(),
                             selected = NULL)
      })
      
      output$debug <- renderPrint(
        values$comb()
      )
    }
  )
}