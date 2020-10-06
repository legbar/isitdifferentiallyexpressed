#Make tidy data.frame
make_tidy <- function(x, first_column){
  as.data.frame(x) %>%
    rownames_to_column(first_column)
}
#summarise res
make_res_summary <- function(res, fdr) {
  res %>%
    filter(padj < fdr) %>%
    mutate(Sign = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated")) %>%
    group_by(Sign) %>%
    summarise(Number = n())
}