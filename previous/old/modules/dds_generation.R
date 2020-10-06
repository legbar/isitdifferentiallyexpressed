#load tidyverse
library(tidyverse)
#Storing gene metadata in the rowData of a SummarizedExperiment: https://support.bioconductor.org/p/62374/
#Kallisto alignments were done with Ensembl v98. TRAP: A cDNA/ncRNA fasta with hSNCA added.
#Salmon quant with Ensembl v101 cDNA fasta only with SNCA cDNA added

#Sample metadata
meta_trap <- readRDS("data/meta_trap.rds")
#Quantification files
files <- file.path("/zfs/analysis/trap/active/alignment/salmon/cohorts_combined/", meta_trap$code, "quant.sf")
names(files) <- meta_trap$name
#tx2gene for tximport
library(ensembldb)
ensDb_mouse <- ensDbFromGtf("ftp://ftp.ensembl.org/pub/release-101/gtf/mus_musculus/Mus_musculus.GRCm38.101.gtf.gz")
ensDb_mouse <- EnsDb(ensDb_mouse)
k_mouse <- keys(ensDb_mouse, keytype = "TXNAME")

tx2gene_mouse <- select(ensDb_mouse, 
                  keys = k_mouse, 
                  columns = c("GENEID"), 
                  keytype = "TXNAME")
tx2gene_mouse <- tx2gene_mouse[,-3]
ensDb_human <- ensDbFromGtf("ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz")
ensDb_human <- EnsDb(ensDb_human)
k_human <- keys(ensDb_human, keytype = "TXNAME")
tx2gene_human <- select(ensDb_human, 
                        keys = k_human, 
                        columns = c("GENEID"), 
                        keytype = "TXNAME")
tx2gene_human <- tx2gene_human[,-3]
tx2gene_SNCA <- tx2gene_human[tx2gene_human$GENEID == "ENSG00000145335", -3]
tx2gene_trap <- rbind(tx2gene_mouse, tx2gene_SNCA)
saveRDS(tx2gene_trap, "data/tx2gene_trap.rds")
#tximport
library(tximport)
txi <- tximport(files, 
                type = "salmon", 
                tx2gene = tx2gene_trap, 
                ignoreTxVersion = T)
#DDS object, removing 0 count genes
library(DESeq2)
dds <- DESeqDataSetFromTximport(txi, 
                                meta_trap, 
                                design = ~1)
saveRDS(dds, "data/dds_trap.rds")

#rowData for dds
library(biomaRt)
ensembl_mouse <- useMart("ensembl", 
                        dataset = "mmusculus_gene_ensembl", 
                        host = "www.ensembl.org")
anno_mouse <- getBM(attributes = c('ensembl_gene_id', 
                                   'external_gene_name', 
                                   'gene_biotype',
                                   'description'), 
                        filters = 'ensembl_gene_id', 
                        values = rownames(dds), 
                        mart = ensembl_mouse)
#human homologs: Some genes have multiple homologs
homolog_mouse <- getBM(attributes = c('ensembl_gene_id', 
                                      'hsapiens_homolog_ensembl_gene',
                                      'hsapiens_homolog_associated_gene_name'), 
                       filters = 'ensembl_gene_id', 
                       values = anno_mouse$ensembl_gene_id, 
                       mart = ensembl_mouse)
homolog_mouse_multiple <- homolog_mouse %>%
  group_by(ensembl_gene_id) %>%
  summarise(n = n()) %>%
  filter(n > 1)
anno_mouse <- left_join(anno_mouse,
                        homolog_mouse, 
                        by = "ensembl_gene_id")
ensembl_human <- useMart("ensembl", 
                         dataset = "hsapiens_gene_ensembl", 
                         host = "www.ensembl.org")
anno_SNCA <- getBM(attributes = c('ensembl_gene_id', 
                                   'external_gene_name', 
                                   'gene_biotype',
                                   'description'), 
                    filters = 'ensembl_gene_id', 
                    values = rownames(dds[c("ENSG00000145335"),]), 
                    mart = ensembl_human) %>%
  mutate(hsapiens_homolog_ensembl_gene = "ENSG00000145335", 
         hsapiens_homolog_associated_gene_name = "SNCA")
anno_trap <- rbind(anno_mouse,
                  anno_SNCA)
saveRDS(anno_trap, "data/anno_trap.rds")


dds <- dds[,(meta_trap$cohort == "PK1" & meta_trap$ip == "ip" & meta_trap$outlier == FALSE)]
dds <- DESeq(dds)
dds <- dds[rowData(dds)$allZero == FALSE]
dim(dds)

saveRDS(dds, "data/dds_PK1_IP.rds")

  

