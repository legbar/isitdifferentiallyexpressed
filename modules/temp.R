






counts_rank <- counts(dds_kw_trap[,colnames(dds_kw_trap) %in% (sample_metadata %>% filter(ip == ip))$sample_name]) %>%
  make_tidy("gene") %>%
  pivot_longer(-gene, 
               names_to = "sample_name", 
               values_to = "count") %>%
  group_by(gene) %>%
  summarise(count = mean(count)) %>%
  arrange(desc(count)) %>%
  mutate(rank = row_number(), 
         trap_enrichment = ifelse(gene %in% dopaminergic_genes$gene, "Enriched", 
                                  ifelse(gene %in% non_dopaminergic_genes$gene, "Depleted", 
                                         "Unchanged"))) %>%
  inner_join(anno_mmus, by = c("gene" = "ensembl_gene_id"))





  select("Gene name" = external_gene_name, 
         "Human Homolog" = hsapiens_homolog_associated_gene_name, 
         "TRAP Status" = trap_enrichment,
         "Rank Expression" = rank, 
         "Description" = description) 



ggplot(counts_rank, 
       aes(x = rank, 
           y = count)) +
  geom_point(aes(colour = trap_enrichment),
             alpha = 1) +
  geom_label_repel(data = subset(counts_rank, external_gene_name %in% c("Th", 
                                                                        "Slc6a3", 
                                                                        "Gad2", 
                                                                        "Actb", 
                                                                        "Lrrk2", 
                                                                        "Gba", 
                                                                        "Aldh1l1", 
                                                                        "Nfe2l2", 
                                                                        "Slc18a2", 
                                                                        "Orai1", 
                                                                        "Stim1", 
                                                                        "Syt1", 
                                                                        "Fasn", 
                                                                        "Snca", 
                                                                        "Gapdh", 
                                                                        "Tfeb", 
                                                                        "Tfe3", 
                                                                        "Rest", 
                                                                        "Nedd4")), 
                   aes(x = rank, 
                       y = count, 
                       label = external_gene_name, 
                       colour = trap_enrichment), 
                   box.padding = 0.5, 
                   force = 3, 
                   key_glyph = "point") +
  geom_hline(yintercept = 10, 
             linetype = "dotted") +
  geom_text(aes(x = 1, 
                y = 10, 
                label = "Filtered", 
                vjust = -0.5)) +
  scale_y_log10(breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000)) +
  labs(y = "Log10 Counts", 
       x = "Ranked Expression", 
       colour = "TRAP Enrichment") +
  theme_cowplot(cowplot_size) +
  scale_x_continuous(trans = "reverse")
# scale_colour_manual(values = c("Red", "Darkgreen", "Lightblue"))
