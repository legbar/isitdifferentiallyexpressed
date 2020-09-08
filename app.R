library(shiny)
library(tidyverse)
library(shinydashboard)
library(data.table)
library(tximport)
library(DESeq2)
library(IHW)
library(plotly)
library(cowplot)
library(ggsci)
library(ggrepel)
library(ggbeeswarm)


meta_trap <- readRDS("data/meta_trap.rds")
anno_mmus <- readRDS("data/anno_mmus.rds")
counts_rank <- readRDS("data/counts_rank.rds")
vsd_kw_all_ip <- readRDS("data/vsd_kw_all_ip.rds")
vsd_pk1_ip_mb <- readRDS("data/vsd_pk1_ip_mb.rds")
age_genes_common <- readRDS("data/age_genes_common.rds")
res_kw_age_ip_anno <- readRDS("data/res_kw_age_ip_anno.rds")
res_pk1_age_ip_anno <- readRDS("data/res_pk1_age_ip_anno.rds")
dds_kw_age_ip <- readRDS("data/dds_kw_age_ip.rds")
dds_pk1_age_ip <- readRDS("data/dds_pk1_age_ip.rds")

setDT(counts_rank)
setDT(meta_trap)

tx2gene <- readRDS("data/tx2gene_trap.rds")

fdr < 0.05

#Custom functions
#Make tidy data.frame
make_tidy <- function(x, first_column){
  as.data.frame(x) %>%
    rownames_to_column(first_column)
}
#summarise res
make_res_summary <- function (res, fdr) {
  res %>%
    filter(padj < fdr) %>%
    mutate(Sign = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated")) %>%
    group_by(Sign) %>%
    summarise(Number = n())
}

header <- dashboardHeader(title = "isitdifferentiallyexpressed?")

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("How expressed is my gene?", tabName = "expression_check", icon = icon("check-circle")),
    menuItem("What does my gene do in ageing?", tabName = "ageing_check", icon = icon("check-circle")),
    menuItem("Custom Analysis", tabName = "metadata_build", icon = icon("table"))
  )
)
 
 body <- dashboardBody(
  tabItems(
    tabItem(tabName = "metadata_build",
            fluidRow(
              column(width = 3, 
                     box(title = "Introduction", 
                         height = NULL, 
                         collapsible = FALSE, 
                         solidHeader = TRUE, 
                         status = "info", 
                         width = NULL, 
                         p("This is an introduction to the website.")
                         ), 
                     selectInput(inputId = "experiment_choice", 
                                 label = "Choose Experiment", 
                                 choices = c("TRAP", "iPSC"), 
                                 selected = NULL,
                                 selectize = TRUE), 
                     conditionalPanel(
                       condition = 'input.experiment_choice == "TRAP"', 
                       selectizeInput(inputId = "trap_choice", 
                                     label = "Which TRAP Cohort?", 
                                     choices = unique(meta_trap$cohort))),
                     conditionalPanel(
                       condition = 'input.experiment_choice == "TRAP"', 
                       checkboxGroupInput(inputId = "age_choice", 
                                      label = "Which age groups?", 
                                      choices = unique(meta_trap$age), 
                                      select = unique(meta_trap$age))),
                     conditionalPanel(
                       condition = 'input.trap_choice != "PK1"', 
                       checkboxGroupInput(inputId = "ovx_choice", 
                                      label = "Which genotypes?", 
                                      choices = unique(meta_trap$ovx), 
                                      select = unique(meta_trap$ovx))),
                     conditionalPanel(
                       condition = 'input.experiment_choice == "TRAP" & input.trap_choice == "PK1"', 
                       checkboxGroupInput(inputId = "region_choice", 
                                      label = "Which regions?", 
                                      choices = unique(meta_trap$region), 
                                      select = unique(meta_trap$region))),
                     conditionalPanel(
                       condition = 'input.experiment_choice == "TRAP"', 
                       checkboxInput(inputId = "ip_only",
                                     label = "TRAP Samples Only", 
                                     value = TRUE)),
                     conditionalPanel(
                       condition = 'input.trap_choice == "PK1"', 
                       checkboxInput(inputId = "remove_outliers", 
                                     label = "Remove outliers", 
                                     value = TRUE)),
                     selectizeInput(inputId = "exclude_samples_choice", 
                                    label = "Exclude any samples?", 
                                    choices = NULL, 
                                    multiple = TRUE), 
                     selectizeInput(inputId = "comparison_choice", 
                                    label = "Comparison Factor", 
                                    choices = NULL),
                     fluidRow(
                       column(
                         width = 4, 
                         selectizeInput(inputId = "nominator_choice", 
                                        label = "Nominator", 
                                        choices = NULL)
                       ), 
                       column(
                         width = 1, 
                         hr(),
                         tags$b("vs")
                       ),
                       column(
                         width = 4, 
                         selectizeInput(inputId = "denominator_choice", 
                                        label = "Denominator", 
                                        choices = NULL)
                       )
                     ),
                     selectizeInput(inputId = "covariates_choice", 
                                    label = "Covariates to include:", 
                                    multiple = TRUE,
                                    choices = NULL),
                     verbatimTextOutput("final_rows_value"),
                     uiOutput("ui_makedds"), 
                     hr(),
                     verbatimTextOutput("debugdiy"), 
                     box(
                       title = "Result Summary", 
                       DT::dataTableOutput("res_summary"), 
                       width = NULL, 
                       solidHeader = TRUE
                     )
                     ),
              column(
                width = 9, 
                box(
                  title = "TRAP Metadata",
                  DT::dataTableOutput("meta_select"), 
                  width = NULL, 
                  solidHeader = TRUE
                  ), 
                box(title = "Results Plot", 
                    plotlyOutput("results_plot"))
                )
              )
            ),
    
    
    #NEXT PAGE
    tabItem(
      tabName = "expression_check", 
      fluidRow(
        column(width = 2, 
               box(title = "Information", 
                   height = NULL, 
                   collapsible = FALSE, 
                   solidHeader = TRUE, 
                   status = "info", 
                   width = NULL, 
                   p("To find the level of expression of a gene of interest, type its name in the box below and click it or hit enter. To remove genes, click it and hit backspace.")), 
               selectizeInput(inputId = "expression_choices",
                              label = "Check Genes", 
                              selected = c("Th", "Gfap"),
                              multiple = TRUE, 
                              choices = counts_rank$external_gene_name)
      ), 
      column(width = 6, 
             box(title = "TRAP Expression Ranked", 
                 plotOutput("counts_rank_plot"), 
                 width = NULL, 
                 solidHeader = TRUE), 
             box(title = "Gene Information", 
                 DT::dataTableOutput("counts_rank_table"), 
                 width = NULL, 
                 solidHeader = TRUE)), 
      column(width = 4, 
             box(title = "TRAP Status Distribution", 
                 plotOutput("counts_rank_density"), 
                 width = NULL, 
                 solidHeader = TRUE), 
             box(title = "Description", 
                 height = NULL, 
                 collapsible = FALSE, 
                 solidHeader = FALSE, 
                 status = "info", 
                 width = NULL, 
                 p("The above plot shows the distribution of counts for genes measured in TRAP. They are split into 3 categories based on their abundance in TRAP samples relative to a brain homogenate reference. There is considerable overlap between the categories, except for the most highly expressed genes."))
      )
    )
    ), 
    tabItem(
      tabName = "ageing_check", 
      fluidRow(
        column(width = 2, 
               selectizeInput(inputId = "ageing_gene_choices",
                              label = "Check Genes", 
                              selected = c("Syt1"),
                              multiple = FALSE, 
                              choices = anno_mmus[anno_mmus$ensembl_gene_id %in% age_genes_common, ]$external_gene_name)), 
        column(width = 4,
               box(title = "Cohort 1", 
                   plotOutput("ageing_counts_plot_kw"), 
                   width = NULL, 
                   solidHeader = TRUE)), 
        column(width = 4,
               box(title = "Cohort 2", 
                   plotOutput("ageing_counts_plot_pk1"), 
                   width = NULL, 
                   solidHeader = TRUE))
              
      ), 
      fluidRow(
        column(width = 4, 
               offset = 2, 
               box(title = NULL, 
               plotOutput("ageing_volcano_plot_kw"), 
               width = NULL, 
               solidHeader = TRUE)), 
        column(width = 4,
               box(title = NULL, 
                   plotOutput("ageing_volcano_plot_pk1"), 
                   width = NULL, 
                   solidHeader = TRUE))
      )
    )
  )
  )

ui <- dashboardPage(header, 
                    sidebar, 
                    body)

server <- function(input, output, session) {
  
  #Make a list of reactive values - not currently in use
  values <- reactiveValues()
  
  #Rows selected by TRAP choice
  trap_choice_rows <- reactive({
    meta_trap[cohort %in% input$trap_choice, which = TRUE]
  })
  #Rows selected by excluding outliers
  outlier_rows <- reactive({
    if (input$remove_outliers == TRUE) {
    meta_trap[outlier == FALSE, which = TRUE]
    } else {
    meta_trap[outlier == TRUE | outlier == FALSE, which = TRUE]
    }
  })
  #Rows selected by choosing only IP samples
  ip_rows <- reactive({
    if (input$ip_only == TRUE) {
      meta_trap[ip == "ip", which = TRUE]
    } else {
      meta_trap[ip %in% meta_trap$ip, which = TRUE]
    }
  })
  #Rows from age_choice
  age_choice_rows <- reactive({
    meta_trap[age %in% input$age_choice, which = TRUE]
  })
  #Rows from OVX choice
  ovx_choice_rows <- reactive({
    meta_trap[ovx %in% input$ovx_choice, which = TRUE]
  })
  #Rows from region choice
  region_choice_rows <- reactive({
    meta_trap[region %in% input$region_choice, which = TRUE]
  })
  
  meta0 <- reactive({
    rows <- Reduce(intersect, 
                   list(trap_choice_rows(),
                        outlier_rows(),
                        ip_rows(), 
                        age_choice_rows(), 
                        ovx_choice_rows(), 
                        region_choice_rows()))
    meta_trap[rows]
  })
  observeEvent(meta0(), {
    updateSelectizeInput(session,
                         "exclude_samples_choice", 
                         choices = meta0()$name, 
                         selected = lapply(reactiveValuesToList(input), unclass)$exclude_samples_choice)
  })
  
  #Rows excluded specifically
  exclude_samples_rows <- reactive({
    meta_trap[!name %in% input$exclude_samples_choice, which = TRUE ]
  })
  
  
  #Build the filtered table by intersecting a list of all the row numbers
  meta <- reactive({
    final_rows <- Reduce(intersect,
                         list(trap_choice_rows(),
                              outlier_rows(),
                              ip_rows(), 
                              age_choice_rows(), 
                              ovx_choice_rows(), 
                              region_choice_rows(), 
                              exclude_samples_rows()))
    meta_trap[final_rows]
  })
  
  #Provide comparison options and create a filter based on whether there is more than one factor left in the filtered metadata
  comparison_options <- c("age", "ovx", "region", "sex")
  comparison_options_filter <- reactive({
    lapply(meta()[, ..comparison_options], function(x) length(unique(x))) > 1
  })
  #Provide covariate options and create a filter based on whether there is more than one factor left in the filtered metadata
  covariates_options <- c("th_enrichment", "age", "sex", "ovx", "region")
  covariates_options_filter <- reactive({
    lapply(meta()[, ..covariates_options], function(x) length(unique(x))) > 1
  })

  comparison_choices_listen <- reactive({
    list(input$trap_choice,
         input$remove_outliers,
         input$age_choice,
         input$ovx_choice, 
         input$region_choice, 
         input$exclude_samples_choice)
  })

  #Update the comparison choices available based on the TRAP experiment selected
  observeEvent(comparison_choices_listen(), {
    updateSelectizeInput(session, 
                         "comparison_choice", 
                         choices = {
                           options <- comparison_options[comparison_options_filter()]
                           }, 
                         selected = lapply(reactiveValuesToList(input), unclass)$comparison_choice)
  })
  
  nominator_choices_listen <- reactive({
    list(input$trap_choice,
         input$remove_outliers,
         input$age_choice,
         input$ovx_choice, 
         input$region_choice, 
         input$exclude_samples_choice, 
         input$comparison_choice)
  })
  
  #Update the nominator choices available based on the comparison selected
  observeEvent(nominator_choices_listen(), {
    updateSelectizeInput(session,
                         "nominator_choice",
                         choices = unique(meta()[[input$comparison_choice]]), 
    selected = lapply(reactiveValuesToList(input), unclass)$nominator_choice)
    })
  
  #Update the demoninator choices available based on the nominator selected
  observeEvent(input$nominator_choice, {
    updateSelectizeInput(session, 
                         "denominator_choice", 
                         choices = {
                           options <- unique(meta()[[input$comparison_choice]])
                           options[!options == input$nominator_choice]
                         }, 
                         selected = lapply(reactiveValuesToList(input), unclass)$denominator_choice)
  })
  

  
  #This listening function allow the covariates options to update when any of the below listed inputs changes. 
  # I could add the other options (age, region), but it might be better to write a function that checks for 
  # errors at the end instead and stops the button actioning.
  covariates_choices_listen <- reactive({
    list(input$trap_choice,
         input$remove_outliers,
         input$age_choice,
         input$ovx_choice, 
         input$region_choice, 
         input$exclude_samples_choice,
         input$comparison_choice)
  })
  
  observeEvent(covariates_choices_listen(), {
    updateSelectizeInput(session, 
                         "covariates_choice", 
                         choices = {
                           options <- covariates_options[covariates_options_filter()]
                           options[options != input$comparison_choice] 
                         },
                         selected = "th_enrichment")
  })
  
  # Filter data based on selections
  output$meta_select <- DT::renderDataTable({
    meta()
  })

  output$ui_makedds <- renderUI({
    actionButton("button_makedds", "Generate the dataset", class = "btn btn-success")
  })
  
  generate_dds <- reactive({
    #baseMean for filtering before independent filtering
    #This influences WGCNA results, because of the resulting matrix size
    basemean <- 10
    
    # files <- file.path("/home/ubuntu/isitdifferentiallyexpressed/input_data/kallisto", meta()$code, "abundance.h5")
    files <- file.path("/zfs/analysis/thesis-pk/input_data/kallisto_PK1_KW2", meta()$code, "abundance.h5")
    names(files) <- meta()$name
    txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = T)

    dds <- DESeqDataSetFromTximport(txi, colData = meta(), design = as.formula(paste0("~", paste(input$covariates_choice, input$comparison_choice, sep = " + "), collapse = " + ")))
    dds[[input$comparison_choice]] <- relevel(dds[[input$comparison_choice]], input$denominator_choice)
    dds <- estimateSizeFactors(dds)
    keep_feature <- rowMeans(counts(dds, normalized = TRUE)) >= basemean
    dds_removed <- dds[!keep_feature, ]
    dds <- dds[keep_feature, ]
    dds <- DESeq(dds,
                 minReplicatesForReplace = Inf,
                 parallel = TRUE)
    
    print(resultsNames(dds))

    return(dds)
  })
  
  generate_counts <- reactive({
    
  })
  
  generate_res <- reactive({
    # Set FDR
    fdr = 0.05
    # Set lfc at 5%
    lfc = log2(1.05)
    #Filtering with independent hypothesis weighting
    filterfunction <- "ihw"
    #Ensuring independent filtering is performed
    independentfiltering = TRUE
    #Using apeglm shrinkage (note this requires use of the coef argument as opposed to contrast for interaction results (see ?results))
    shrinkage = "apeglm"
    
    res <- results(values$dds,
                   contrast = c(input$comparison_choice, input$nominator_choice, input$denominator_choice),
                   alpha = fdr,
                   lfcThreshold = lfc,
                   filterFun = get(filterfunction),
                   independentFiltering = independentfiltering,
                   parallel = TRUE)
    
    res_shrunken <- lfcShrink(values$dds,
                              coef = paste(input$comparison_choice,
                                           input$nominator_choice,
                                           "vs",
                                           input$denominator_choice, sep = "_"),
                              type = "apeglm",
                              lfcThreshold = lfc,
                              res = res, quiet = T)

    res$log2FoldChange <- res_shrunken$log2FoldChange
    
    res <- make_tidy(res, "gene") %>% 
      mutate(significant = ifelse(padj < fdr, TRUE, FALSE), 
             padj = -log10(padj)) %>% 
      inner_join(anno_mmus, by = c("gene" = "ensembl_gene_id"))
      
    return(res)
  })
  
  generate_res_summary <- reactive({
    # Set FDR
    res_summary <- make_res_summary(values$res, 0.05)
    return(res_summary)
  })
  
  observeEvent(input$button_makedds,
               {
                 values$dds <- generate_dds()
                 values$res <- generate_res()
                 values$res_summary <- generate_res_summary()
               })
  
  output$res_summary <- DT::renderDataTable({
    if(!is.null(values$res_summary)){
      values$res_summary
    }
  })
  
  # output$results_plot <- renderPlotly({
  #   if(!is.null(values$res_summary)){
  #     g <- make_tidy(values$res, "gene") %>%
  #     ggplot(aes(y = -log10(padj), 
  #                x = log2FoldChange, 
  #                label = gene)) +
  #     geom_point(alpha = 0.25, 
  #                size = 2) +
  #       theme_cowplot() +
  #       scale_color_lancet()
  #     # return(g)
  #     ggplotly(g)
  #   }
  # })
  
  output$results_plot <- renderPlotly({
    if(!is.null(values$res_summary)){
      # plot_ly(data = values$res) %>%
      #   add_trace(y = ~padj, 
      #             x = ~log2FoldChange,
      #             type = "scattergl", 
      #             mode = "markers",
      #             marker = list(size = 12), 
      #             hoverinfo = "text", 
      #             text = paste("Gene Name: ", values$res[["external_gene_name"]], "\n", 
      #                          "Gene Type: ", values$res[["gene_biotype"]], 
      #                          sep = ""),
      #             alpha = 0.3) %>%
      #   layout(showlegend = FALSE, 
      #          )
      
      p <- ggplot(values$res, 
             aes(x = log2FoldChange, 
                 y = padj, 
                 colour = significant, 
                 label = external_gene_name, 
                 text = paste("Gene Name: ", external_gene_name, sep = ""))) +
        geom_point(alpha = 0.25, 
                   size = 2) +
        theme_cowplot(10) +
        scale_color_d3() +
        labs(x = "Log2 fold change",
             y = "-Log<sub>10</sub> P adjusted value") +
        theme(legend.position = "none")
      fig <- ggplotly(p, tooltip = "text") %>%
        toWebGL()
      fig
      }
  })
  
  
  output$counts_rank_plot <- renderPlot({
    p <- ggplot(counts_rank, 
             aes(x = rank, 
                 y = count)) +
      geom_point(aes(colour = trap_enrichment),
                 alpha = 0.1, 
                 position = position_jitter(h = 0.1, w = 0.1), size = 1) +
      geom_hline(yintercept = 10, 
                 linetype = "dotted") +
      geom_text(aes(x = 1, 
                    y = 10, 
                    label = "Filtered", 
                    vjust = -0.5), 
                size = 6) +
      geom_label_repel(data = subset(counts_rank, external_gene_name %in% input$expression_choices), 
                       aes(x = rank, 
                           y = count, 
                           label = external_gene_name, 
                           colour = trap_enrichment), 
                       box.padding = 0.5, 
                       force = 3, 
                       key_glyph = "point", 
                       fontface = "bold", 
                       size = 8) +
      scale_y_log10(breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000)) +
      labs(y = expression("Log"[10]~"Counts"), 
           x = "Ranked Expression", 
           colour = "TRAP Status") +
      theme_cowplot() +
      scale_x_continuous(trans = "reverse", limits = c(NA, 1)) +
      scale_color_manual(values = c("Green", "Orange", "Red")) +
      theme(legend.position = "top", 
            legend.justification = c(0, 1))
    
    # fig <- ggplotly(p) %>%
    #   toWebGL()
    # fig
    
    p
    
  })
  
  output$counts_rank_density <- renderPlot({
    p <- ggplot(counts_rank, 
                aes(x = log10(count + 1), 
                fill = trap_enrichment, 
                colour = trap_enrichment)) +
      geom_density(alpha = 0.5) +
      geom_vline(data = subset(counts_rank, external_gene_name %in% input$expression_choices), 
                 aes(xintercept = log10(count + 1), 
                     colour = trap_enrichment), 
                 linetype = "dotted") +
      geom_text(data = subset(counts_rank, external_gene_name %in% input$expression_choices), 
                aes(x = log10(count + 1), 
                    y = Inf, 
                    label = external_gene_name,
                    colour = trap_enrichment), 
                vjust = 1, 
                hjust = -0.1, 
                size = 3)  +
      scale_color_manual(values = c("Green", "Orange", "Red")) +
      scale_fill_manual(values = c("Green", "Orange", "Red")) +
      theme_cowplot() +
      theme(legend.position = "top") +
      labs(colour = "", 
           fill = "", 
           x = expression("Log"[10]~"Counts"), 
           y = "Density")
      
    p
  })
  

  output$counts_rank_table <- DT::renderDataTable({
    counts_rank %>%
      filter(external_gene_name %in% input$expression_choices) %>%
      select("Gene name" = external_gene_name, 
             "Human Homolog" = hsapiens_homolog_associated_gene_name, 
             "TRAP Status" = trap_enrichment,
             "Rank Expression" = rank, 
             "Description" = description) 
    })
  
  ageing_counts_plot_function <- function(dds){
    counts(dds, normalized = TRUE) %>%
    make_tidy("gene") %>%
    pivot_longer(-gene, names_to = "name", values_to = "count") %>%
    inner_join(meta_trap, by = "name") %>%
    inner_join(anno_mmus, by = c("gene" = "ensembl_gene_id")) %>%
      mutate(age = recode(age, old = "Old", young = "Young")) %>%
    filter(external_gene_name %in% input$ageing_gene_choices) %>%
    ggplot(aes(x = age, 
               y = count, 
               # colour = age, 
               fill = age)) +
      geom_boxplot() +
    geom_point(position = position_jitterdodge()) +
    expand_limits(x = 0, y = 0) +
    theme_cowplot() +
    scale_color_d3() +
      scale_fill_d3() +
      labs(x = "Age", 
           y = "Count", 
           colour = "Age")
  }
  
  output$ageing_counts_plot_kw <- renderPlot({
    if(!is.null(input$ageing_gene_choices)){
    ageing_counts_plot_function(dds_kw_age_ip)
    }
  })
  
  output$ageing_counts_plot_pk1 <- renderPlot({
    if(!is.null(input$ageing_gene_choices)){
    ageing_counts_plot_function(dds_pk1_age_ip)
    }
  })
  
  ageing_volcano_plot_function <- function(res_anno){
    fdr < 0.05
    p <- res_anno %>%
      ggplot(aes(label = external_gene_name, 
                 colour = ifelse(padj < fdr, TRUE, FALSE))) +
      geom_point(aes(x = log2FoldChange, 
                     y = -log10(padj)), 
                 alpha = 0.25, 
                 size = 2) +
      geom_point(data = subset(res_anno, external_gene_name %in% input$ageing_gene_choices),
                 aes(x = log2FoldChange, 
                     y = -log10(padj)), 
                 size = 3) +
      geom_label_repel(data = subset(res_anno, external_gene_name %in% input$ageing_gene_choices),
                       aes(x = log2FoldChange, 
                           y = -log10(padj))) +
      theme_cowplot() +
      scale_color_d3() +
      labs(x = "Log2 fold change",
           y = "-Log<sub>10</sub> P adjusted value") +
      theme(legend.position = "none")
    
    p
  }
  
  output$ageing_volcano_plot_kw <- renderPlot({
    ageing_volcano_plot_function(res_kw_age_ip_anno)
  })
  output$ageing_volcano_plot_pk1 <- renderPlot({
    ageing_volcano_plot_function(res_pk1_age_ip_anno)
  })
  
    output$results_plot <- renderPlot({
         ageing_volcano_plot_function(res_kw_age_ip) 
    })
  

    output$debugdiy <- renderPrint({
    if(!is.null(values$res)){
      print(values$dds)
      print(resultsNames(values$dds))
      print(input$comparison_choice)
      print(input$nominator_choice)
      print(input$denominator_choice)
      print(summary(values$res))
    }
  })
  
  
  # output$dds_design <- renderUI({
  #   if(is.null(values$exp_design))
  #     return(NULL)
  #   covars <- colnames(values$exp_design)
  #   selectInput("dds_design", 
  #               label = "Select covariates to model: ", 
  #               choices = c(NULL, covars), 
  #               selected = NULL, 
  #               multiple = TRUE)
  # })
  
}

shinyApp(ui, server)