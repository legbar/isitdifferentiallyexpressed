library(shiny)
library(tidyverse)
library(shinydashboard)
# library(data.table)
# library(plotly)
# library(cowplot)
# library(ggsci)
# library(ggrepel)

# dds_trap <- readRDS("data/dds_trap.rds")
# anno_trap <- readRDS("data/anno_trap.rds")

header <- dashboardHeader(title = "isitdifferentiallyexpressed?",
                          titleWidth = 300)

sidebar <- dashboardSidebar(width = 300,
                            sidebarMenu(
                              # menuItem("How expressed is my gene?", tabName = "expression_check", icon = icon("check-circle")),
                              # menuItem("What does my gene do in ageing?", tabName = "ageing_check", icon = icon("check-circle")),
                              menuItem(
                                "Custom Dataset",
                                tabName = "customDataTab",
                                icon = icon("table")
                              ))
)


body <- dashboardBody(tabItems(tabItem(
  tabName = "customDataTab",
  tabsetPanel(
    id = "customDataTabset",
    tabPanel("Build Dataset",
             fluidRow(
               column(
                 width = 4,
                 box(
                   title = "Configure Dataset",
                   c(
                     # customDataBuildSpecificUI("customDataBuildSpecific"),
                     customDataBuildGenericUI("customDataBuild")
                   ),
                   width = TRUE,
                   solidHeader = TRUE
                 )
               ),
               column(
                 width = 8,
                 box(
                   title = "Sample Information",
                   DT::dataTableOutput("metaFiltered"),
                   width = NULL,
                   solidHeader = TRUE
                 )
               )
             )),
    tabPanel("View Results",
             customDataViewUI("customDataView"),)
  )
)))

ui <- dashboardPage(header,
                    sidebar,
                    body)


server <- function(input, output, session) {
  # customDataBuildSpecificValues <-
  #   customDataBuildSpecificServer("customDataBuildSpecific")
  # customDataBuildGenericValues <-
  #   customDataBuildGenericServer("customDataBuildGeneric", customDataBuildSpecificValues)
  
  customDataBuildValues <- customDataBuildGenericServer("customDataBuild")
  
  # customDataViewValues <-
  #   customDataViewServer("customDataView", customDataBuildValues)
  
  
  # Filter data based on selections
  output$metaFiltered <- DT::renderDataTable({
    customDataBuildValues$metaFiltered() %>%
      select("Name" = name,
             "Cohort" = cohort,
             Region,
             Age,
             Genotype)
  })
}


shinyApp(ui, server)

