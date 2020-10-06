library(shiny)
library(shinyjs)
library(tidyverse)
library(shinydashboard)
library(DT)
library(data.table)
library(tximport)
library(DESeq2)
library(cowplot)
library(ggsci)
library(ggrepel)
library(plotly)

header <- dashboardHeader(title = "Is it differentially expressed?",
                          titleWidth = 300)

sidebar <- dashboardSidebar(width = 175,
                            sidebarMenu(
                              menuItem(
                                "Custom Dataset",
                                tabName = "customDataTab",
                                icon = icon("table")
                              )
                            ))

body <- dashboardBody(
  useShinyjs(),
  tags$head(
    tags$link(rel = "stylesheet", 
              type = "text/css", 
              href = "custom.css")),
  tabItems(tabItem(
  tabName = "customDataTab",
  tabsetPanel(id = "customDataTabset",
              tabPanel(title = "Build Dataset",
                        value = "build",
                       fluidRow(
                         column(
                           width = 3,
                           box(
                             title = "Configure Dataset",
                             ddsBuild_skeleton_UI("ddsBuild_skeleton"),
                             width = TRUE,
                             solidHeader = TRUE
                           )
                         ), 
                         column(
                           width = 9, 
                           box(
                             title = "Sample Information", 
                             DT::dataTableOutput("metaFinal"), 
                             width = TRUE, 
                             solidHeader = TRUE
                           ), 
                           box(
                             title = "DEBUG", 
                             verbatimTextOutput("debug"),
                             width = TRUE, 
                             solidHeader = TRUE
                           )
                         )
                       )) 
              )
)))


ui <- dashboardPage(header,
                    sidebar,
                    body)


server <- function(input, output, session) {
  
  addClass(selector = "body", class = "sidebar-collapse")
  
  hideTab(inputId = "customDataTabset", 
          target = "viewRes")
  
  values <- ddsBuild_skeleton_Server("ddsBuild_skeleton")
  
  output$metaFinal <- DT::renderDataTable({
    values$metaFinal()[,values$displayColumns()]
  })
  
  observeEvent(values$dds(), {

    removeTab(inputId = "customDataTabset", 
              target = "resBuildView")
    
    appendTab(inputId = "customDataTabset",
              tab = tabPanel(title = "View results", 
                        value = "resBuildView",
                        fluidRow(
                          resBuildView_UI("resBuildView")
                        )))
  })
  
  test <- resBuildView_Server("resBuildView", values)

    output$debug <- renderPrint({
      req(values$dds())
      test$comparisonType()
    })
}


shinyApp(ui, server)
