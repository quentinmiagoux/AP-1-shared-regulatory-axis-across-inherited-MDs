# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with DESeq2
library(DESeq2)

# Set to TRUE only when intentionally regenerating data/EDMD_paper.RDS.
write_processed_data <- FALSE

# Load counts table and annotations from GEO.
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE204804", "file=GSE204804_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")

# load gene annotations 
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID

# sample selection
gsms <- "111111000000XXXXXX000XXX000XXXXXX000"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
tbl <- tbl[ ,sel]

# group membership for samples
gs <- factor(sml)
groups <- make.names(c("LMNA","WT"))
levels(gs) <- groups
sample_info <- data.frame(Group = gs, row.names = colnames(tbl))

# pre-filter low count genes
# keep genes with at least N counts > 10, where N = size of smallest group
keep <- rowSums( tbl >= 10 ) >= min(table(gs))
tbl <- tbl[keep, ]

ds <- DESeqDataSetFromMatrix(countData=tbl, colData=sample_info, design= ~Group)

ds <- DESeq(ds, test="Wald", sfType="poscount")

# extract results for top genes table
r <- results(ds, contrast=c("Group", groups[1], groups[2]), alpha=0.05, pAdjustMethod ="fdr")

tT <- r[order(r$padj),] 
tT <- merge(as.data.frame(tT), annot, by=0, sort=F)

tT <- subset(tT, select=c("GeneID","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean","Symbol","Description"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

################################################################
#   General expression data visualization
dat <- log2(counts(ds, normalized = T) + 1) # extract normalized counts


ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)

test <- getBM(attributes = c('gene_biotype','entrezgene_id','hgnc_symbol'),
              mart = ensembl) %>% 
  dplyr::filter(gene_biotype == "protein_coding") %>% 
  dplyr::select(-gene_biotype)


data <- tT %>%
  left_join(test,by=c("Symbol"="hgnc_symbol")) %>% 
  filter(!is.na(entrezgene_id)) %>% arrange(desc(log2FoldChange))


data_filter <- data %>% filter(padj < 0.05, abs(log2FoldChange) > log2(1.5)) 
data_h_degs_up <- data_filter %>% filter(log2FoldChange > 0)
data_h_degs_down <- data_filter %>% filter(log2FoldChange < 0)

gsea_all_genes_h <- data$log2FoldChange
names(gsea_all_genes_h) <- data$entrezgene_id
gsea_all_genes_h <- gsea_all_genes_h[!is.na(gsea_all_genes_h)]

# source("Enrichment.R")
# 
# EDMD <- enrich_it_baby(gsea_full =  gsea_all_genes_h, 
#                         deg_down = data_h_degs_down$entrezgene_id, 
#                         deg_up = data_h_degs_up$entrezgene_id, 
#                         universe = data$entrezgene_id)
# 

dat <- dat %>% as.data.frame %>%  rownames_to_column(var='entrezid') %>% mutate(entrezid = as.integer(entrezid))

data_sym <- dat %>%
  left_join(test, by = c("entrezid" = "entrezgene_id")) %>%
  filter(!is.na(hgnc_symbol)) %>% dplyr::select(-entrezid)

# Agrégation par symbole (moyenne)
data_sym_agg <- data_sym %>%
  group_by(hgnc_symbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "hgnc_symbol")




EDMD <- list(dataset = data, expr = data_sym_agg, group = gs)


if (write_processed_data) {
  saveRDS(EDMD, file = "data/EDMD_paper.RDS")
}

