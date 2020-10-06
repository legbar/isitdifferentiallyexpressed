# TRAP DDS Build Questions
ddsBuild_iPSC_UI <- function(id) {
  ns <- NS(id)
  
  meta_trap <- readRDS("data/meta_trap.rds")
  # tx2gene_trap <- readRDS("data/tx2gene_trap.rds")
  
  tagList(
    selectizeInput(
      inputId = ns("ipscDatasetSelect"),
      label = "Which iPSC Dataset?",
      choices = unique(meta_trap$cohort)
    )
  )
}


#Server functions
ddsBuild_iPSC_Server <- function(id) {
  moduleServer(id,
               ## Below is the module function
               function(input, output, session) {
                 
                 
               }
               )
}


  
  
  
