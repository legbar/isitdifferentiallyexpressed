# Res Build/View
resBuildView_UI <- function(id) {
  ns <- NS(id)
  tagList(
    column(
      width = 4,
      box(
        title = "Customise Results",
        tagList(
          selectizeInput(
            inputId = ns("contrastSelect"),
            label = "Select a specific comparison",
            choices = NULL,
            selected = NULL
          )
        ),
        width = TRUE,
        solidHeader = TRUE
      )
    ),
  column(
      width = 8,
      box(
        title = "Results Table",
        DT::dataTableOutput(ns("resultsTable")),
        width = TRUE,
        solidHeader = TRUE
      ),
      
      box(
        title = "Volcano",
        plotlyOutput(ns("volcanoPlot")),
        width = TRUE,
        solidHeader = TRUE
      ),
      
      box(
        title = "Bar",
        plotlyOutput(ns("barPlot")),
        width = TRUE,
        solidHeader = TRUE
      ),
      
      tabBox(
        title = "",
        width = NULL,
        tabPanel(title = tagList(icon("bar-chart"), "Count Distribution"),
                 uiOutput(ns(
                   "sampleDistributionBoxPanel"
                 ))),
        tabPanel(title = tagList(icon("area-chart"), "Density Plot"),
                 uiOutput(ns(
                   "sampleDistributionDensityPanel"
                 ))),
        tabPanel(title = tagList(icon("object-group"), "PCA"),
                 uiOutput(ns("pcaUI")))
      )
  )
  )
}
#Server functions
resBuildView_Server <- function(id, inputValues) {
  moduleServer(id,
               function(input, output, session) {
                 
                 values <- inputValues

                 values$contrastOptions <- reactive({
                   req(values$dds())
                   generateComparisonCombinations(values)
                   })
                 
                 observe({
                   updateSelectizeInput(session,
                                        "contrastSelect",
                                        choices = names(values$contrastOptions()))
                 })
                 
                 values$contrastSelect <- reactive(values$contrastOptions()[input$contrastSelect])
                 
                 values$res <- reactive(
                                             {
                                               if (values$comparisonType() == "FACTOR"){
                                                 results(values$dds(),
                                                         unname(values$contrastSelect())[[1]],
                                                         # lfcThreshold,
                                                         # alpha,
                                                         tidy = TRUE)       
                                               }
                                               }
                                             )
                 
                 output$resultsTable <- DT::renderDataTable({
                   DT::datatable({
                     values$res() %>%
                     select("Gene" = row, baseMean, padj, log2FoldChange)
                   }) %>%
                     formatRound(c("padj", "baseMean", "log2FoldChange"), 3)
                 })
                 
                 output$sampleDistributionBoxPanel <- renderUI({
                   ns <- session$ns
                   tagList(
                     fluidRow(
                       column(
                         width = 3, 
                         selectizeInput(inputId = ns("colourBoxplot"), 
                                        label = "Colour by...", 
                                        choices = values$comparisonOptions(),
                                        selected = values$comparisonSelect()
                         )
                       ), 
                       column(
                         width = 9, 
                         plotOutput(ns("sampleDistributionBox"))
                       )
                     )
                   )
                 })
                 
                 output$sampleDistributionBox <- renderPlot({
                   values$vsd_long_meta() %>%
                     ggplot(aes_string(x = "name", 
                                y = "count", 
                                fill = names(values$comparisonOptions()[values$comparisonOptions() == input$colourBoxplot])), 
                            alpha = 0.5) +
                     geom_boxplot() +
                     theme_cowplot() +
                     theme(axis.text.x = element_text(angle = 90)) +
                     labs(y = expression("Log"[2]~"Counts"), 
                          x = "Sample") +
                     scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16))
                 })
                 
                 output$sampleDistributionDensityPanel <- renderUI({
                   ns <- session$ns
                   tagList(
                     fluidRow(
                       column(
                         width = 3, 
                         selectizeInput(inputId = ns("colourDensity"), 
                                        label = "Colour by...", 
                                        choices = values$comparisonOptions(),
                                        selected = values$comparisonSelect()
                         )
                       ), 
                       column(
                         width = 9, 
                         plotOutput(ns("sampleDistributionDensity"))
                       )
                     )
                   )
                 })
                 
                 output$sampleDistributionDensity <- renderPlot({
                   values$vsd_long_meta() %>%
                     ggplot(aes_string(x = "count", 
                                       fill = "name",
                                       colour = names(values$comparisonOptions()[values$comparisonOptions() == input$colourDensity])), 
                            alpha = 0.5) +
                     geom_line(stat = "density") +
                     theme_cowplot() +
                     theme(axis.text.x = element_text(angle = 90)) +
                     labs(y = expression("Log"[2]~"Counts"), 
                          x = "Sample")
                 })

                 output$pcaUI <- renderUI({
                   ns <- session$ns
                   tagList(
                     fluidRow(
                       column(
                         width = 3, 
                         selectizeInput(inputId = ns("colourPCA"), 
                                        label = "Colour by...", 
                                        choices = values$comparisonOptions(),
                                        selected = values$comparisonSelect()
                         ),
                         selectizeInput(inputId = ns("shapePCA"), 
                                        label = "Shape by...", 
                                        choices = values$comparisonOptions(),
                                        selected = values$comparisonSelect()
                         )
                       ), 
                       column(
                         width = 9, 
                         plotOutput(ns("pcaBiplot"))
                       )
                     )
                   )
                 })
                 
                 output$pcaBiplot <- renderPlot({
                   
                   values$pca_meta() %>%
                     ggplot(aes_string(x="PC1",
                                       y="PC2", 
                                       # label = "name", 
                                       color = names(values$comparisonOptions()[values$comparisonOptions() == input$colourPCA]), 
                                       shape = names(values$comparisonOptions()[values$comparisonOptions() == input$shapePCA]))) +
                     geom_point(size = 4) +
                     labs(x=paste0("PC1: ",round(values$pca_percentVar()[1]*100,1),"%"),
                          y=paste0("PC2: ",round(values$pca_percentVar()[2]*100,1),"%")) +
                     # geom_label_repel(box.padding = 0.5) +
                     theme_cowplot() +
                     theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank()) 
                   # ggtitle(title) +
                   # theme(legend.position = "none") +
                   # coord_fixed()
                   
                 })

                 output$volcanoPlot <- renderPlotly({
                   
                   plot_ly(
                     data = values$res(),
                     x = ~ log2FoldChange,
                     y = ~ -log10(values$res()[["padj"]]),
                     type = "scatter",
                     mode = "markers",
                     # color = ~ x,
                     # colors = c(input$downColor, "black", input$upColor),
                     # marker = list(size = input$volcanoPointSize),
                     # hoverinfo = "text",
                     # text = ~ paste(
                     #   "</br>Gene:",
                     #   resultTable()$gene_id,
                     #   "</br>A value:",
                     #   round(a.value, 4),
                     #   "</br>M value:",
                     #   round(m.value, 4),
                     #   "</br>p-value:",
                     #   round(p.value, 4),
                     #   "</br>q-value:",
                     #   round(q.value, 4),
                     #   "</br>Rank:",
                     #   rank
                     # ),
                     key =  ~ values$res()[["row"]],
                     source = "volcano"
                   )
                   
                 })
                 
                 output$barPlot <- renderPlotly({
                   # Read in hover data
                   eventdata <- event_data("plotly_hover", source = "volcano")
                   
                   validate(need(
                     !is.null(eventdata),
                     "Hover over the point to show gene's expression level of interest."
                   ))
                   # Get point number
                   gene_id <- eventdata$key
                   #Get expression values
                   counts <- values$counts() %>%
                     filter(gene == gene_id)
                     
                   
                   # data <- variables$CountData
                   # data.cl <- variables$groupListConvert
                   # 
                   # expression <- t(expression[data.cl != 0])
                   # data.cl <- data.cl[data.cl != 0]
                   # 
                   # xOrder <-
                   #   data.frame("name" = row.names(expression), "group" = data.cl)
                   # xOrderVector <- unique(xOrder[order(xOrder$group),]$name)
                   # xform <- list(categoryorder = "array",
                   #               categoryarray = xOrderVector,
                   #               title = "")
                   
                   plot_ly(
                     x = ~ counts$name,
                     y = ~ counts$count,
                   #   color = as.factor(data.cl),
                   #   text = expression[, 1],
                   #   textposition = 'outside',
                   #   showlegend = FALSE,
                     type = "bar",
                   #   name = "Raw"
                   # ) %>%
                   #   # add_trace(
                   #   #   y = ~ expressionNor[, 1],
                   #   #   text = round(expressionNor[, 1], 2),
                   #   #   textposition = 'auto',
                   #   #   name = "Normalized",
                   #   #   type = "scatter",
                   #   #   mode = "lines+markers",
                   #   #   marker = list(
                   #   #     size = 10,
                   #   #     line = list(color = 'rgba(0, 0, 0, 0)',
                   #   #                 width = 2)
                   # #   )
                   # # ) %>%
                   # layout(
                   #   xaxis = xform,
                   #   yaxis = list(title = "Raw Count"),
                   #   title = colnames(expression)
                   )
                 })
                 
                 
                 return(values)
  })
}
