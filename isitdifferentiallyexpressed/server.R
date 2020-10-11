

dds_snca <-
    readRDS("R/dds_snca.rds")
dds_lrrk2 <-
    readRDS(
        "R/dds_lrrk2.rds"
    )
dds_gba <-
    readRDS("R/dds_gba.rds")
anno_human <-
    readRDS(
        "R/anno_human.rds"
    )

results_list <- readRDS(
    "R/results_list.rds"
)
comparison_choices <- readRDS(
    "R/comparison_choices.rds"
)
comparison_dictionary <- readRDS(
    "R/comparison_dictionary.rds"
)

signif_check <- readRDS(
    "R/signif_check.rds"
)



# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    observeEvent(input$dds_select, {
        updateSelectizeInput(session,
                             inputId = "comparison_select",
                             choices = comparison_choices[[input$dds_select]])
    })
    
    dds_select <- reactive({
        eval(parse(text = input$dds_select))
    })
    
    results <- reactive({
        req(input$comparison_select %in% comparison_choices[[input$dds_select]])

        results_list[[input$dds_select]][[input$comparison_select]]
        
    })
    
    results_filtered <- reactive({
        req(results())
        
        if (input$signif_only == TRUE){
            rows <- which(results()$padj < 0.1)
        } else {
            rows <- seq(1, nrow(results()), by = 1)
        }
        
        results()[rows,] %>%
            filter(baseMean > input$baseMean_select)
    })
    
    output$results_table <- DT::renderDataTable({
        req(input$comparison_select %in% comparison_choices[[input$dds_select]])
        results_filtered() %>%
            mutate(
                Direction = ifelse(
                    sign(log2FoldChange) == 1,
                    "Upregulated",
                    "Downregulated"
                ),
                # Significant = ifelse(padj < input$alpha_select, "Significant", "NS"),
                padj = signif(padj, 3),
                log2FoldChange = signif(log2FoldChange, 3),
                baseMean = signif(baseMean, 5)
            ) %>%
            select(
                "Gene" = external_gene_name,
                "Log 2 Fold Change" = log2FoldChange,
                "Adjusted P value" = padj,
                "Mean counts" = baseMean,
                "Description" = description
            )
    },
    rownames = FALSE,
    selection = "single")
    
    gene_select <- reactive({
        req(input$results_table_rows_selected != "")
        results_filtered()[input$results_table_rows_selected, ]$ensembl_gene_id
    })

    output$scatterPlot <- renderPlot({
        req(input$comparison_select %in% comparison_choices[[input$dds_select]])
        
            g <- ggplot({
                counts(dds_select(), normalized = TRUE)[rownames(dds_select()) == gene_select(), ] %>%
                    as_tibble(rownames = "Name") %>%
                    inner_join(as_tibble(colData(dds_select())), by = "Name")
            },
            aes(
                x = eval(parse(text = comparison_dictionary[[input$dds_select]][["contrasts"]][[input$comparison_select]][1])),
                y = value,
                colour = eval(parse(text = comparison_dictionary[[input$dds_select]][["contrasts"]][[input$comparison_select]][1]))
            )) +
                theme_cowplot() +
                labs(colour = comparison_dictionary[[input$dds_select]][["contrasts"]][[input$comparison_select]][1],
                     x = comparison_dictionary[[input$dds_select]][["contrasts"]][[input$comparison_select]][1],
                     y = "Count") +
                scale_y_continuous(limits = c(1, NA))
            if (str_detect(input$comparison_select, "vs")) {
                g + 
                    geom_point() +
                    geom_boxplot() 
        } else {
            g +
                geom_point()
        }
        
    })
    
    # output$debug <- ({
    #     renderPrint({
    #     })
    # })
    
    updateSelectizeInput(session, 
                         'byGene_gene_select', 
                         choices = unique(signif_check$external_gene_name), 
                         selected = "SNCA",
                         server = TRUE)
    
    output$signif_check <- renderText({
        data <- signif_check %>%
            filter(external_gene_name == input$byGene_gene_select)
        
        paste("<p><b>", data$dataset_comparison, "</b>", ": ", ifelse(data$significant == TRUE, "<font color=\"#FF0000\"><b> Differentially expressed </b></font>", 
                                                                   ifelse(data$significant == FALSE, "Not differentially expressed", 
                                                                          "Not quantified due to low counts in some samples")), "</p>", sep = "")
        
    })
    
    updateSelectizeInput(session, 
                         'byGene_comparison_select', 
                         choices = unique(signif_check$dataset_comparison), 
                         selected = NULL,
                         server = TRUE)

    byGene_reactives <- reactive({
        {
            req(input$byGene_gene_select, input$byGene_comparison_select)
            if (str_detect(input$byGene_comparison_select, "SNCA")) {
                return(list("dds_snca", 
                            counts(dds_snca, normalized = TRUE)[anno_human[which(anno_human$external_gene_name == input$byGene_gene_select),]$ensembl_gene_id,], 
                            results_list[["dds_snca"]][[sub(".*? ", "", input$byGene_comparison_select)]] %>% filter(external_gene_name == input$byGene_gene_select)))
            }
            if (str_detect(input$byGene_comparison_select, "LRRK2")) {
                return(list("dds_lrrk2", 
                            counts(dds_lrrk2, normalized = TRUE)[anno_human[which(anno_human$external_gene_name == input$byGene_gene_select),]$ensembl_gene_id,], 
                            results_list[["dds_lrrk2"]][[sub(".*? ", "", input$byGene_comparison_select)]] %>% filter(external_gene_name == input$byGene_gene_select)))
            }
            if (str_detect(input$byGene_comparison_select, "GBA")) {
                return(list("dds_gba", 
                            counts(dds_gba, normalized = TRUE)[anno_human[which(anno_human$external_gene_name == input$byGene_gene_select),]$ensembl_gene_id,], 
                            results_list[["dds_gba"]][[sub(".*? ", "", input$byGene_comparison_select)]] %>% filter(external_gene_name == input$byGene_gene_select)))
            }
        }
    })
    
    output$byGene_scatterPlot <- renderPlot({
        req(byGene_reactives())

        g_byGene <- ggplot({
            byGene_reactives()[[2]] %>%
                as_tibble(rownames = "Name") %>%
                inner_join(as_tibble(colData(eval(parse(text = byGene_reactives()[[1]])))), by = "Name")
        },
        aes(
            x = eval(parse(text = comparison_dictionary[[byGene_reactives()[[1]]]][["byGene_meta_column"]][[input$byGene_comparison_select]])),
            y = value,
            colour = eval(parse(text = comparison_dictionary[[byGene_reactives()[[1]]]][["byGene_meta_column"]][[input$byGene_comparison_select]])))) +
            theme_cowplot() +
            labs(colour = comparison_dictionary[[byGene_reactives()[[1]]]][["byGene_meta_column"]][[input$byGene_comparison_select]],
                 x = comparison_dictionary[[byGene_reactives()[[1]]]][["byGene_meta_column"]][[input$byGene_comparison_select]],
                 y = "Count") +
            scale_y_continuous(limits = c(1, NA))

        if (str_detect(input$byGene_comparison_select, "vs")) {
            g_byGene + 
                geom_boxplot() +
                geom_point()
        } else {
            g_byGene +
                geom_point()
        }

    })
    
    # output$byGene_debug <- 
    #     renderPrint({
    #        byGene_reactives()[[3]]
    #         })
    
    output$byGene_resultTable <-
        DT::renderDataTable({
            datatable({byGene_reactives()[[3]] %>%
                    mutate(log2FoldChange = signif(log2FoldChange, 3), 
                           padj = signif(padj, 3), 
                           baseMean = signif(baseMean, 3)) %>%
                    select("Gene" = external_gene_name, 
                           "Log 2 Fold Change" = log2FoldChange, 
                           "Adjusted P" = padj, 
                           "Mean counts" = baseMean)}, options = list(dom = 't'), rownames = FALSE)
        })

    waiter_hide()
})
