

generateDDS <-
  function(files,
           meta,
           TX2GENE,
           quantType = "kallisto",
           covariates,
           comparison,
           basemean = 10) {
    incProgress(0.25, detail = paste("Gathering samples"))
    
    names(files) <- meta$name
    txi <-
      tximport(
        files,
        type = quantType,
        tx2gene = TX2GENE,
        ignoreTxVersion = T
      )
    
    dds <-
      DESeqDataSetFromTximport(txi,
                               colData = meta,
                               design = as.formula(paste0(
                                 "~",
                                 paste(covariates, comparison, sep = " + "),
                                 collapse = " + "
                               )))
    
    incProgress(0.25, detail = paste("Normalising data"))
    
    # dds[[input$comparisonSelect]] <- relevel(dds[[input$comparisonSelect]], input$denominatorSelect)
    dds <- estimateSizeFactors(dds)
    keep_feature <-
      rowMeans(counts(dds, normalized = TRUE)) >= basemean
    dds_removed <- dds[!keep_feature, ]
    dds <- dds[keep_feature, ]
    dds <- DESeq(dds,
                 minReplicatesForReplace = Inf,
                 parallel = TRUE)
    
    incProgress(0.25, detail = paste("Filtering low count genes"))
    Sys.sleep(1)
    
    incProgress(0.25, detail = paste("Ready for viewing"))
    Sys.sleep(2)
    
    print("DDS Object build complete.")
    
    return(dds)
    
  }

generateComparisonCombinations <- function(values) {
  combinations <- combn(unique(values$metaFinal()[[values$comparisonSelect()]]), 2, simplify = F)
  combinations_names <- combn(unique(values$metaFinal()[[names(values$comparisonOptions()[values$comparisonOptions() == values$comparisonSelect()])]]), 2, simplify = F)
  contrasts <- c(lapply(combinations, function(x){c(as.character(values$comparisonSelect()), as.character(x[1]), as.character(x[2]))}),
                 lapply(combinations, function(x){c(as.character(values$comparisonSelect()), as.character(x[2]), as.character(x[1]))})) 

  contrasts_names <- c(lapply(combinations_names, function(x){paste(names(values$comparisonOptions()[values$comparisonOptions() == values$comparisonSelect()]), paste(x, collapse = " versus "), sep = ": ")}),
                       lapply(combinations_names, function(x){paste(names(values$comparisonOptions()[values$comparisonOptions() == values$comparisonSelect()]), paste(rev(x), collapse = " versus "), sep = ": ")})) 

  names(contrasts) <- contrasts_names

  return(contrasts)
  
  }
  
generateRES_factor <- function() {
  
  # dds, contrast, lfcThreshold, alpha
  results(values$dds(),
          input$contrastSelect,
          # lfcThreshold,
          # alpha,
          tidy = TRUE)
}