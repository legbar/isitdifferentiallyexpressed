

dds_snca <-
    readRDS("/zfs/analysis/isitdifferentiallyexpressed_objectCreation/output/dds_snca.rds")
dds_lrrk2 <-
    readRDS(
        "/zfs/analysis/isitdifferentiallyexpressed_objectCreation/output/dds_lrrk2.rds"
    )
dds_gba <-
    readRDS("/zfs/analysis/isitdifferentiallyexpressed_objectCreation/output/dds_gba.rds")

dds_TRAP_ALL <-
    readRDS("/zfs/analysis/isitdifferentiallyexpressed_objectCreation/output/dds_TRAP_ALL.rds")
dds_TRAP_IP <-
    readRDS("/zfs/analysis/isitdifferentiallyexpressed_objectCreation/output/dds_TRAP_IP.rds")

anno_human <-
    readRDS(
        "/zfs/analysis/isitdifferentiallyexpressed_objectCreation/output/anno_human.rds"
    )

anno_mouse <-
  readRDS(
    "/zfs/analysis/isitdifferentiallyexpressed_objectCreation/output/anno_mouse.rds"
  )

results_list <- readRDS(
    "/zfs/analysis/isitdifferentiallyexpressed_objectCreation/output/results_list.rds"
)

comparison_dictionary <- readRDS(
    "/zfs/analysis/isitdifferentiallyexpressed_objectCreation/output/comparison_dictionary.rds"
)


signif_check_human <- readRDS(
    "/zfs/analysis/isitdifferentiallyexpressed_objectCreation/output/signif_check_human.rds"
)
signif_check_mouse <- readRDS(
    "/zfs/analysis/isitdifferentiallyexpressed_objectCreation/output/signif_check_mouse.rds"
)
# will need to import these objects into local directory!
query_genes <- readRDS(
    "/zfs/analysis/isitdifferentiallyexpressed_objectCreation/output/query_genes.rds"
)

human_to_mouse <- readRDS(
    "/zfs/analysis/isitdifferentiallyexpressed_objectCreation/output/human_to_mouse.rds"
) %>%
    distinct() %>%
    rename("homolog" = "mmusculus_homolog_associated_gene_name") %>%
  add_row(external_gene_name = "S100B", 
          homolog = "S100b")

mouse_to_human <- readRDS(
    "/zfs/analysis/isitdifferentiallyexpressed_objectCreation/output/mouse_to_human.rds"
) %>%
    distinct()  %>%
    rename("homolog" = "hsapiens_homolog_associated_gene_name") %>%
  add_row(external_gene_name = "S100b", 
          homolog = "S100B")

enrich_summary <- readRDS(
    "/zfs/analysis/isitdifferentiallyexpressed_objectCreation/output/enrich_summary.rds"
)

comparison_to_dds <- readRDS(
  "/zfs/analysis/isitdifferentiallyexpressed_objectCreation/output/comparison_to_dds.rds"
)

# Define server logic 
shinyServer(function(input, output, session) {
    
    # per dataset tab
    # observe dataset selection and update comparison selection based on that
    observeEvent(input$dds_select, {
        updateSelectizeInput(session,
                             inputId = "comparison_select",
                             choices = names(comparison_dictionary[[input$dds_select]][["contrasts"]]))
    })
    
    # define dds_select as the dds object selected
    dds_select <- reactive({
        eval(parse(text = input$dds_select))
    })
    
    # define results as the result of dds-comparison choice
    results <- reactive({
        req(input$comparison_select %in% names(comparison_dictionary[[input$dds_select]][["contrasts"]]))

        results_list[[input$dds_select]][[input$comparison_select]]
        
    })
    
    # filter results table by either signif, all or baseMean
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
    
    # render results table 
    output$results_table <- DT::renderDataTable({
        req(input$comparison_select %in% names(comparison_dictionary[[input$dds_select]][["contrasts"]]))
        results_filtered() %>%
            mutate(
                # Significant = ifelse(padj < input$alpha_select, "Significant", "NS"),
                log2FoldChange = signif(log2FoldChange, 3), 
                lfcSE = signif(lfcSE, 3),
                pvalue = signif(pvalue, 3),
                padj = signif(padj, 3), 
                baseMean = signif(baseMean, 5)
            ) %>%
            select(
                "Gene" = external_gene_name, 
                "Log<sub>2</sub> Fold Change" = log2FoldChange, 
                "Log<sub>2</sub> Fold Change Standard Error" = lfcSE,
                "Raw P" = pvalue,
                "Adjusted P" = padj, 
                "Mean counts" = baseMean,
                "Description" = description
            )
    },
    rownames = FALSE,
    escape = FALSE,
    selection = "single")
    
    # define gene_select, based on clicking table row
    gene_select <- reactive({
        req(input$results_table_rows_selected != "")
        results_filtered()[input$results_table_rows_selected, ]$ensembl_gene_id
    })

    byDataset_counts <- reactive({
      counts(dds_select(), normalized = TRUE)[rownames(dds_select()) == gene_select(), ] %>%
        as_tibble(rownames = "Name") %>%
        inner_join(colData(dds_select()), copy = TRUE)
    })
    
    
    # render a scatterplot for gene_select
    # obtain normalized counts for the gene of interest in this dds object
    # add colData to counts
    # set x axis to the first term in the result contrast (this wont work for interaction results)
    # set y to count value
    # set colour to the first term in the result contrast
    
    # this is where counts(dds) and colData(dds) is needed
    output$scatterPlot <- renderPlot({
        req(input$comparison_select %in% names(comparison_dictionary[[input$dds_select]][["contrasts"]]))

            g <- ggplot(
              byDataset_counts(),
              aes_string(
                x = comparison_dictionary[[input$dds_select]][["byDataset_meta_column"]][[input$comparison_select]],
                y = "value",
                colour = comparison_dictionary[[input$dds_select]][["byDataset_meta_column"]][[input$comparison_select]])) +
                theme_cowplot() +
                labs(colour = comparison_dictionary[[input$dds_select]][["byDataset_meta_column"]][[input$comparison_select]],
                     x = comparison_dictionary[[input$dds_select]][["byDataset_meta_column"]][[input$comparison_select]],
                     y = "Count") 
                # scale_y_continuous(limits = c(1, NA))
            if (str_detect(input$comparison_select, "vs") | str_detect(input$dds_select, "TRAP")) {
                g +
                    geom_jitter(width = 0.2)
                    # geom_boxplot()
        } else {
            g +
                geom_point()
        }

    })
    
    
    
    # step 1. a gene (human or mouse) gets chosen
    updateSelectizeInput(session, 
                         'byGene_gene_select', 
                         # choices = unique(signif_check$external_gene_name), 
                         choices = query_genes,
                         selected = "SNCA",
                         server = TRUE)
    
    # step 2. a check is made to see what species the gene is from
    human_gene <- reactive({
      req(input$byGene_gene_select)
      
        if (input$byGene_gene_select %in% signif_check_human$external_gene_name) {
            return(TRUE)
        } else {
            return(FALSE)
        }
    })
    
 

    # step 3. the homolog(s) of the alternate species are retrieved
    homologs <- reactive({
      req(input$byGene_gene_select)
      
        if (human_gene() == TRUE) {
          return(human_to_mouse[human_to_mouse$external_gene_name == input$byGene_gene_select,])
          } else {
          return(mouse_to_human[mouse_to_human$external_gene_name == input$byGene_gene_select,])
          }
      })
    
    
  
    # # step 4. if there is more than one homolog, an option is given to select one
    output$human_homolog_select <- renderUI({
      req(homologs())
        if (!human_gene()) {
          if (nrow(homologs()) > 1){
            selectizeInput("human_homolog_select", "Choose a homolog", homologs()$homolog)
          } else if (nrow(homologs()) == 0){
            paste("No homolog exists or was tested")
          } else {
            paste("Homolog: ", homologs()$homolog)
          }
        }
    })

    output$byGene_debug <-
      renderPrint({
        homologs()
      })
    
    output$mouse_homolog_select <- renderUI({
      req(homologs())
      if (human_gene()) {
        if (nrow(homologs()) > 1){
          selectizeInput("mouse_homolog_select", "Choose a homolog", homologs()$homolog)
        } else if (nrow(homologs()) == 0){
          paste("No homolog exists or was tested")
        } else {
          paste("Homolog: ", homologs()$homolog)
        }
      }
    })

    # # step 5. the active gene relationship is established
    active_human_gene <- reactive({
      req(homologs())
        if (human_gene()){
            return(input$byGene_gene_select)
        } else if (nrow(homologs()) > 1) {
            return(input$human_homolog_select)
        } else if (nrow(homologs()) == 1) {
            return(homologs()$homolog)
        } else {
            return(NA)
        }
    })

    active_mouse_gene <- reactive({
      req(homologs())

        if (!human_gene()){
            return(input$byGene_gene_select)
        } else if (nrow(homologs()) > 1) {
            return(input$mouse_homolog_select)
        } else if (nrow(homologs()) == 1) {
            return(homologs()$homolog)
        } else {
            return(NA)
        }
    })

    # filter the comparisons
    human_comparisons <- reactive({
      req(active_human_gene())

        if (is.na(active_human_gene())) {
            return(NULL)
        } else {
            signif_check_human %>%
                filter(external_gene_name == active_human_gene()) %>%
                return()
        }

    })

    mouse_comparisons <- reactive({
      req(active_mouse_gene())

        if (is.na(active_mouse_gene())) {
            return(NULL)
        } else {
            signif_check_mouse %>%
                filter(external_gene_name == active_mouse_gene()) %>%
                return()
        }

    })

    mouse_enrichment <- reactive({
        req(active_mouse_gene())
          if (is.na(active_mouse_gene())) {
            return(NULL)
        } else {
            enrich_summary %>%
                filter(external_gene_name == active_mouse_gene()) %>%
                return()
        }

    })

    output$signif_check_human <- renderText({
      req(human_comparisons())

      paste("<p><b>", human_comparisons()$dataset_comparison, "</b>", ": ", ifelse(human_comparisons()$significant == TRUE, "<font color=\"#1C8E3E\"><b> Differentially expressed </b></font>",
                                                                      ifelse(human_comparisons()$significant == FALSE, "Not differentially expressed",
                                                                             "Not quantified due to low counts in some samples")), "</p>", sep = "")

    })



    output$enrichment_check_mouse <- renderText({
      req(mouse_enrichment())
        paste("<p><b>", mouse_enrichment()$dataset_comparison, "</b>", ": ", ifelse(mouse_enrichment()$enrichment == "enriched", "<font color=\"#1C8E3E\"><b> Enriched </b></font>",
                                                                ifelse(mouse_enrichment()$enrichment == "depleted", "Depleted",
                                                                       "Not differentially expressed")), "</p>")
    })

    output$signif_check_mouse <- renderText({
      req(mouse_comparisons())
        paste("<p><b>", mouse_comparisons()$dataset_comparison, "</b>", ": ", ifelse(mouse_comparisons()$significant == TRUE, "<font color=\"#1C8E3E\"><b> Differentially expressed </b></font>",
                                                                                     ifelse(mouse_comparisons()$significant == FALSE, "Not differentially expressed",
                                                                                            "Not quantified due to low counts in some samples")), "</p>", sep = "")

    })


    # remove per gene plot choices if a homolog for that gene is unavailable
    perGene_plot_choices <- reactive({
      req(active_mouse_gene(), active_human_gene())
        if (is.na(active_mouse_gene())){
            return(unique(signif_check_human$dataset_comparison))
        } else if (is.na(active_human_gene())){
            return(unique(c(enrich_summary$dataset_comparison,
                   signif_check_mouse$dataset_comparison)))
        } else {
            return(unique(c(signif_check_human$dataset_comparison,
                            enrich_summary$dataset_comparison,
                     signif_check_mouse$dataset_comparison)))}})

observe({
  
  updateSelectizeInput(session,
                       'byGene_comparison_select',
                       choices = perGene_plot_choices(),
                       # choices = comparison_to_dds$name,
                       selected = NULL,
                       server = TRUE)
})
  


 





    byGene_dds <- reactive({
      eval(parse(text = comparison_to_dds[comparison_to_dds$name == input$byGene_comparison_select,]$dds))
    })



    selected_gene <- reactive({
      req(active_mouse_gene(), active_human_gene())
      if (str_detect(input$byGene_comparison_select, "TRAP")){
        selected_gene <- anno_mouse[anno_mouse$external_gene_name == active_mouse_gene(),]$ensembl_gene_id
      } else {
        selected_gene <- anno_human[anno_human$external_gene_name == active_human_gene(),]$ensembl_gene_id
      }
      
      
    })
      
    byGene_counts <- reactive({
      req(byGene_dds())
      req(active_mouse_gene(), active_human_gene())
      

  

      counts(byGene_dds(), normalized = TRUE)[rownames(byGene_dds()) == selected_gene(), ] %>%
        as_tibble(rownames = "Name") %>%
        inner_join(colData(byGene_dds()), copy = TRUE)
    })

    # this is where counts(dds) and colData(dds) is needed
    output$byGene_scatterPlot <- renderPlot({
      req(byGene_counts())

      g <- ggplot(
        byGene_counts(),
        aes_string(
          x = comparison_dictionary[[comparison_to_dds[comparison_to_dds$name == input$byGene_comparison_select,]$dds]][["byGene_meta_column"]][[input$byGene_comparison_select]],
          y = "value",
          colour = comparison_dictionary[[comparison_to_dds[comparison_to_dds$name == input$byGene_comparison_select,]$dds]][["byGene_meta_column"]][[input$byGene_comparison_select]])) +

    theme_cowplot() +
    labs(colour = comparison_dictionary[[comparison_to_dds[comparison_to_dds$name == input$byGene_comparison_select,]$dds]][["byDataset_meta_column"]][[input$byGene_comparison_select]],
         x = comparison_dictionary[[comparison_to_dds[comparison_to_dds$name == input$byGene_comparison_select,]$dds]][["byDataset_meta_column"]][[input$byGene_comparison_select]],
         y = "Count")
    # scale_y_continuous(limits = c(1, NA))
      if (!str_detect(input$byGene_comparison_select, "Ageing")) {
        g +
          geom_jitter(width = 0.2)
        # geom_boxplot()
      } else {
        g +
          geom_point()
      }

    })
    
    # output$debug <- ({
    #     renderPrint({
    #       selected_gene()
    #     })
    # })
    
    output$byGene_resultTable <-
      DT::renderDataTable({
      req(byGene_counts())
        datatable({results_list[[comparison_to_dds[comparison_to_dds$name == input$byGene_comparison_select,]$dds]][[sub(".*? ", "", input$byGene_comparison_select)]] %>%
            filter(ensembl_gene_id == selected_gene()) %>%
            mutate(log2FoldChange = signif(log2FoldChange, 3),
                   lfcSE = signif(lfcSE, 3),
                   pvalue = signif(pvalue, 3),
                   padj = signif(padj, 3),
                   baseMean = signif(baseMean, 3)) %>%
            select("Gene" = external_gene_name,
                   "Log<sub>2</sub> Fold Change" = log2FoldChange,
                   "Log<sub>2</sub> Fold Change Standard Error" = lfcSE,
                   "Raw P" = pvalue,
                   "Adjusted P" = padj,
                   "Mean counts" = baseMean)}, options = list(dom = 't'), rownames = FALSE, escape = FALSE)
      })
    
# 
#     byGene_reactives <- reactive({
#         {
#             req(input$byGene_gene_select, input$byGene_comparison_select, active_mouse_gene(), active_human_gene())
#             if (str_detect(input$byGene_comparison_select, "SNCA")) {
#                 return(list("dds_snca", 
#                             counts(dds_snca, normalized = TRUE)[anno_human[which(anno_human$external_gene_name == active_human_gene()),]$ensembl_gene_id,], 
#                             results_list[["dds_snca"]][[sub(".*? ", "", input$byGene_comparison_select)]] %>% filter(external_gene_name == active_human_gene())))
#                             # assay(vst(dds_snca))[anno_human[which(anno_human$external_gene_name == input$byGene_gene_select),]$ensembl_gene_id,]))
#             }
#             if (str_detect(input$byGene_comparison_select, "LRRK2")) {
#                 return(list("dds_lrrk2", 
#                             counts(dds_lrrk2, normalized = TRUE)[anno_human[which(anno_human$external_gene_name == active_human_gene()),]$ensembl_gene_id,], 
#                             results_list[["dds_lrrk2"]][[sub(".*? ", "", input$byGene_comparison_select)]] %>% filter(external_gene_name == active_human_gene())))
#                             # assay(vst(dds_lrrk2))[anno_human[which(anno_human$external_gene_name == input$byGene_gene_select),]$ensembl_gene_id,]))
#             }
#             if (str_detect(input$byGene_comparison_select, "GBA")) {
#                 return(list("dds_gba", 
#                             counts(dds_gba, normalized = TRUE)[anno_human[which(anno_human$external_gene_name == active_human_gene()),]$ensembl_gene_id,], 
#                             results_list[["dds_gba"]][[sub(".*? ", "", input$byGene_comparison_select)]] %>% filter(external_gene_name == active_human_gene())))
#                             # assay(vst(dds_lrrk2))[anno_human[which(anno_human$external_gene_name == input$byGene_gene_select),]$ensembl_gene_id,]))
#             }
#             if (str_detect(input$byGene_comparison_select, "TRAP Enrichment")) {
#                 return(list("dds_TRAP_ALL", 
#                             counts(dds_TRAP_ALL, normalized = TRUE)[anno_mouse[which(anno_mouse$external_gene_name == active_mouse_gene()),]$ensembl_gene_id,], 
#                             results_list[["dds_TRAP_ALL"]][[sub(".*? ", "", input$byGene_comparison_select)]] %>% filter(external_gene_name == active_mouse_gene())))
#                             # assay(vst(dds_lrrk2))[anno_human[which(anno_human$external_gene_name == input$byGene_gene_select),]$ensembl_gene_id,]))
#             } else {
#                 return(list("dds_TRAP_IP", 
#                             counts(dds_TRAP_ALL, normalized = TRUE)[anno_mouse[which(anno_mouse$external_gene_name == active_mouse_gene()),]$ensembl_gene_id,], 
#                             results_list[["dds_TRAP_IP"]][[sub(".*? ", "", input$byGene_comparison_select)]] %>% filter(external_gene_name == active_mouse_gene())))
#             }
#             
#         }
#     })
    
    # output$byGene_scatterPlot <- renderPlot({
    #     req(byGene_reactives())
    #     
    #     generate_ggplot <- function(count_data){
    #         g <- as_tibble(count_data, rownames = "Name") %>%
    #             inner_join(as_tibble(colData(eval(parse(text = byGene_reactives()[[1]])))), by = "Name") %>%
    #             ggplot(
    #                 aes(
    #                     x = eval(parse(text = comparison_dictionary[[byGene_reactives()[[1]]]][["byGene_meta_column"]][[input$byGene_comparison_select]])),
    #                     y = value,
    #                     colour = eval(parse(text = comparison_dictionary[[byGene_reactives()[[1]]]][["byGene_meta_column"]][[input$byGene_comparison_select]])))) +
    #             theme_cowplot() +
    #             labs(colour = comparison_dictionary[[byGene_reactives()[[1]]]][["byGene_meta_column"]][[input$byGene_comparison_select]],
    #                  x = comparison_dictionary[[byGene_reactives()[[1]]]][["byGene_meta_column"]][[input$byGene_comparison_select]])
    #         return(g)
    #     }
    #     
    #     if (str_detect(input$byGene_comparison_select, "vs")) {
    #         byGene_reactives()[[2]] %>%
    #          generate_ggplot() + 
    #             # geom_boxplot() +
    #             geom_point() +
    #             labs(y = "Count") +
    #             scale_y_continuous(limits = c(1, NA))
    #     } else {
    #         if(input$brent_mode == TRUE){
    #             byGene_reactives()[[4]] %>%
    #                 generate_ggplot() + 
    #                 geom_point() +
    #                 labs(y = "Log2 Count") +
    #                 # geom_smooth(method = "lm", se = F) +
    #                 stat_smooth(method = "lm", se = F) +
    #                 stat_regline_equation(aes(label = paste(..rr.label..)))
    #         } else {
    #             byGene_reactives()[[2]] %>%
    #                 generate_ggplot() + 
    #                 geom_point() +
    #                 labs(y = "Count")
    #         }
    #         
    #     }
    #     
        
    # 
    # })
  
  
    


    waiter_hide()
})
