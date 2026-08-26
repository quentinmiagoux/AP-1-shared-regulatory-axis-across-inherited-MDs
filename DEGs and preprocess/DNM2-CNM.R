library(GEOquery)
library(tidyverse)
library(DESeq2)
library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(Orthology.eg.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
set.seed(123)

# Set to TRUE only when intentionally regenerating data/DNM2_paper.RDS.
write_processed_data <- FALSE

GSE160078 <- getGEO("GSE160078",GSEMatrix=TRUE)

info <- GSE160078$GSE160078_series_matrix.txt.gz

info_f <- data.frame(condition = info$`treatment:ch1`, sex = info$`Sex:ch1`, tissue = info$`tissue:ch1`, GT = info$`genotype:ch1`) %>% mutate(condition_f = paste0(condition,"_",GT))

raw_counts_path <- "data/GSE160078_Raw_gene_counts_matrix_Cohort_DNM2.txt"
if (!file.exists(raw_counts_path)) {
  stop("Missing raw count matrix: ", raw_counts_path)
}
data <- read.delim(raw_counts_path, check.names = TRUE)


info_f <- info_f %>% filter(condition_f %in% c("ASO Control_Dnm2SL/+", "ASO Control_WT")) %>% dplyr::arrange(desc(condition_f))
data2 <- dplyr::select(data, WT.ASO.Control.1:WT.ASO.Control.4, `Dnm2SL/+.ASO.Control.1`:`Dnm2SL/+.ASO.Control.6`)
rownames(data2) <- data$Ensembl.Gene.ID


# rownames(data) <- data$Ensembl.Gene.ID

dds <- DESeqDataSetFromMatrix(countData = data2,
                              colData = info_f,
                              design= ~ condition_f)
dds <- DESeq(dds)

dat <- log2(counts(dds, normalized = T) + 1)

res <- results(dds, contrast=c('condition_f',unique(info_f$condition_f)[2],unique(info_f$condition_f)[1]), pAdjustMethod ="fdr")


ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)

test <- getBM(attributes = c('ensembl_gene_id', 'mgi_symbol', 'gene_biotype','entrezgene_id'),
              mart = ensembl) %>% dplyr::filter(gene_biotype == "protein_coding") %>% dplyr::select(-gene_biotype)

data <- res %>% data.frame %>% rownames_to_column(var="ensembl_gene_id") %>%
  left_join(test,by=c("ensembl_gene_id"="ensembl_gene_id")) %>%
  dplyr::filter(mgi_symbol != "") %>% 
  group_by(mgi_symbol) %>% arrange(padj) %>% distinct(mgi_symbol, .keep_all = T)


gsea_all_genes <- data$log2FoldChange
names(gsea_all_genes) <- data$ENTREZID

gsea_all_genes <- gsea_all_genes[order(gsea_all_genes, decreasing = T)]


gsea_all_genes <- gsea_all_genes[!is.na(gsea_all_genes)]

library(Orthology.eg.db)

mapfun <- function(mousegenes){
  gns <- mapIds(org.Mm.eg.db, mousegenes, "ENTREZID", "SYMBOL")
  mapped <- biomaRt::select(Orthology.eg.db, gns, "Homo_sapiens","Mus_musculus")
  naind <- is.na(mapped$Homo_sapiens)
  hsymb <- mapIds(org.Hs.eg.db, as.character(mapped$Homo_sapiens[!naind]), "SYMBOL", "ENTREZID")
  out <- data.frame(Mouse_symbol = mousegenes, mapped, Human_symbol = NA)
  out$Human_symbol[!naind] <- hsymb
  out
}


data_annot <- mapfun(data$mgi_symbol)


data_h <- data %>% 
  left_join(., dplyr::select(data_annot, Mouse_symbol,Human_symbol), by =c("mgi_symbol"="Mouse_symbol")) %>% 
  filter(!is.na(Human_symbol)) %>% 
  filter(!is.na(padj))

human_entrezid <-  bitr(data_h$Human_symbol, fromType = "SYMBOL", toType = "ENTREZID",OrgDb = org.Hs.eg.db, drop = F) 

data_h <- left_join(data_h, human_entrezid, by=c("Human_symbol"="SYMBOL")) %>% arrange(desc(log2FoldChange))


gsea_all_genes_h <- data_h$log2FoldChange
names(gsea_all_genes_h) <- data_h$ENTREZID
gsea_all_genes_h <- gsea_all_genes_h[!is.na(gsea_all_genes_h)]




data_h_degs <- data_h %>% filter(padj <= 0.05 & abs(log2FoldChange) >= log2(1.5))

data_h_degs_up <- data_h_degs %>% filter(log2FoldChange > 0)
data_h_degs_down <- data_h_degs %>% filter(log2FoldChange < 0)

# source("Enrichment.R")
# 
# DMN2 <- enrich_it_baby(gsea_full =  gsea_all_genes_h, 
#                deg_down = data_h_degs_down$ENTREZID, 
#                deg_up = data_h_degs_up$ENTREZID, 
#                universe = data_h$ENTREZID)


dat_tmp <- dat %>% data.frame %>% rownames_to_column(var = "ensembl_gene_id") %>%
  left_join(test,by=c("ensembl_gene_id" = "ensembl_gene_id")) %>%
  dplyr::filter(mgi_symbol != "") %>% 
  group_by(mgi_symbol)

data_sym <- dat_tmp %>% 
  left_join(., dplyr::select(data_annot, Mouse_symbol,Human_symbol), by = c("mgi_symbol" = "Mouse_symbol")) %>% 
  filter(!is.na(Human_symbol)) 


data_sym_tmp  <- data_sym %>% 
  ungroup %>% 
  dplyr::select(!c(ensembl_gene_id,mgi_symbol,entrezgene_id)) %>% 
  dplyr::select(Human_symbol, everything())

# Agregation par symbole (moyenne)
data_sym_agg <- data_sym_tmp %>%
  group_by(Human_symbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "Human_symbol")

groups <- data_sym_agg %>% colnames %>% str_split('\\.', simplify = T) %>% .[,1]

DNM2 <- list(dataset = data_h, expr = data_sym_agg, group = groups)

if (write_processed_data) {
  saveRDS(DNM2, file = "data/DNM2_paper.RDS")
}
