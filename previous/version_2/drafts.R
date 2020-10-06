tx2gene_trap <- readRDS("data/tx2gene_trap.rds")
files <- file.path("/zfs/analysis/kallisto_store/", meta_trap$code, "abundance.h5")
names(files) <- meta_trap$name
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene_trap, ignoreTxVersion = T)
dds <- DESeqDataSetFromTximport(txi, 
                                meta_trap,
                                design = ~1)
dds <- DESeq(dds)
counts <- counts(dds, normalized = T) %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "name", values_to = "count") %>%
  inner_join(meta_trap, by = "name") 

counts %>%
  ggplot(aes_string(x = "name", 
                    y = "count", 
                    fill = "Region"), 
         alpha = 0.5) +
  geom_boxplot() +
  scale_y_log10() +
  theme_cowplot() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  xlab("Sample") +
  ylab("Gene Count")

vst(dds, blind = FALSE) %>%
  assay(.) %>%
  as_tibble(rownames = "gene")


list <- list("A" = "a", "B" = "b", "list" = list("first" = "word", "second" = TRUE, "third" = c(1, 2, 3)))
names(list[list == "a"])

list$list[["second"]]








req(input$makeVolcanoPlot)
isolate({
  dt <- resultTable()
  
  downCut <- input$CutFC[1]
  upCut <- input$CutFC[2]
  
  dt$color <- "None"
  tryCatch({
    dt[dt$m.value <= downCut,]$color <- "Down"
    dt[dt$m.value >= upCut,]$color <- "Up"
    dt[dt[[yaxis]] > input$Cutpvalue,]$color <-
      "None"
  }, error = function(e) {
    sendSweetAlert(session = session, title = "ERROR", text = "No data was satisfied to your cut-off!")
  })
  
  x <- factor(dt$color)
  levels(x) <- list("Down" = 0,
                    "None" = 1,
                    "Up" = 2)
  
  # Add annotation
  key <- resultTable()$gene_id
  
  if (is.null(input$resultTableInVolcanalPlot_rows_selected)) {
    annotation <- list()
  } else {
    markerSelect <- dt[input$resultTableInVolcanalPlot_rows_selected, ]
    
    annotation <- list(
      x = markerSelect$m.value,
      y = -log10(markerSelect[[yaxis]]),
      text = markerSelect$gene_id,
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 7,
      ax = 20,
      ay = 40
    )
  }
  
  p <- plot_ly(
    data = dt,
    x = ~ m.value,
    y = ~ -log10(dt[[yaxis]]),
    type = "scatter",
    mode = "markers",
    color = ~ x,
    colors = c(input$downColor, "black", input$upColor),
    marker = list(size = input$volcanoPointSize),
    hoverinfo = "text",
    text = ~ paste(
      "</br>Gene:",
      resultTable()$gene_id,
      "</br>A value:",
      round(a.value, 4),
      "</br>M value:",
      round(m.value, 4),
      "</br>p-value:",
      round(p.value, 4),
      "</br>q-value:",
      round(q.value, 4),
      "</br>Rank:",
      rank
    ),
    key =  ~ key,
    source = "volcano"
  ) %>%
    layout(
      xaxis = list(title = input$xlabs),
      yaxis = list(title = input$ylabs),
      title = input$graphicTitle,
      legend = list(
        orientation = 'h',
        xanchor = "center",
        x = 0.5,
        y = 1.05
      ),
      annotations = annotation,
      shapes = list(
        list(
          type = 'line',
          y0 =  ~ min(-log10(dt[[yaxis]])),
          y1 =  ~ max(-log10(dt[[yaxis]])),
          x0 = upCut,
          x1 = upCut,
          line = list(dash = 'dot', width = 2)
        ),
        list(
          type = 'line',
          y0 =  ~ min(-log10(dt[[yaxis]])),
          y1 =  ~ max(-log10(dt[[yaxis]])),
          x0 = downCut,
          x1 = downCut,
          line = list(dash = 'dot', width = 2)
        ),
        list(
          type = 'line',
          y0 = -log10(input$Cutpvalue),
          y1 = -log10(input$Cutpvalue),
          x0 =  ~ min(m.value),
          x1 =  ~ max(m.value),
          line = list(dash = 'dot', width = 2)
        )
      )
    )
  variables$VolcanoPlotObject <- p
  p
})
})
runVolcano$runVolcanoValue <- input$makeVolcanoPlot
})

# Render volcanoUI ----
output$volcanoUI <- renderUI({
  if(runVolcano$runVolcanoValue){
    tagList(
      fluidRow(
        column(8, plotlyOutput("volcanoPloty") %>% withSpinner()),
        column(4, plotlyOutput("geneBarPlotInVolcano") %>% withSpinner())
      )
    )
  } else {
    helpText("Please click [Generate Volcano Plot] first.")
  }
})



# This function render a button of R code of making vocalno plot ----

observeEvent(input$makeVolcanoPlot, {
  output$runVolcanoPlot <- renderText({
    variables$runVolcanoPlot
  })
})


# This function render a plotly of specific gene expression value in barplot.----

output$geneBarPlotInVolcano <- renderPlotly({
  # Read in hover data
  eventdata <- event_data("plotly_hover", source = "volcano")
  validate(need(
    !is.null(eventdata),
    "Hover over the point to show gene's expression level of interest."
  ))
  # Get point number
  gene_id <- eventdata$key
  # Get expression level (Original)
  expression <-
    variables$CountData[row.names(variables$CountData) == gene_id,]
  # Get expression level (Normalized)
  expressionNor <-
    t(t(variables$norData[row.names(variables$norData) == gene_id,]))
  
  data <- variables$CountData
  data.cl <- variables$groupListConvert
  
  expression <- t(expression[data.cl != 0])
  data.cl <- data.cl[data.cl != 0]
  
  xOrder <-
    data.frame("name" = row.names(expression), "group" = data.cl)
  xOrderVector <- unique(xOrder[order(xOrder$group),]$name)
  xform <- list(categoryorder = "array",
                categoryarray = xOrderVector,
                title = "")
  
  plot_ly(
    x = ~ row.names(expression),
    y = ~ expression[, 1],
    color = as.factor(data.cl),
    text = expression[, 1],
    textposition = 'outside',
    showlegend = FALSE,
    type = "bar",
    name = "Raw"
  ) %>%
    # add_trace(
    #   y = ~ expressionNor[, 1],
    #   text = round(expressionNor[, 1], 2),
    #   textposition = 'auto',
    #   name = "Normalized",
    #   type = "scatter",
    #   mode = "lines+markers",
    #   marker = list(
    #     size = 10,
    #     line = list(color = 'rgba(0, 0, 0, 0)',
    #                 width = 2)
  #   )
  # ) %>%
  layout(
    xaxis = xform,
    yaxis = list(title = "Raw Count"),
    title = colnames(expression)
  )
})
