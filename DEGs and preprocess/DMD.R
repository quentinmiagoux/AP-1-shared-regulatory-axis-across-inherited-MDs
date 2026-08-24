library(tidyverse)
library(biomaRt)
library(clusterProfiler)
library(ReactomePA)

# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE38417", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "1111110000000000000000"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("DMD","WT"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=fit2$genes %>% nrow)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol"))
write.table(tT, file=stdout(), row.names=F, sep="\t")



ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)

test <- getBM(attributes = c('gene_biotype','entrezgene_id','hgnc_symbol'),
              mart = ensembl) %>% dplyr::filter(gene_biotype == "protein_coding") %>% dplyr::select(-gene_biotype)


data <- tT %>% filter( Gene.symbol != "") %>%
  left_join(test,by=c("Gene.symbol"="hgnc_symbol")) %>% 
  filter(!is.na(entrezgene_id)) %>% arrange(desc(logFC))

data_filter <- data %>% filter(adj.P.Val < 0.05, abs(logFC) > log2(1.5))
data_h_degs_up <- data_filter %>% filter(logFC > 0)
data_h_degs_down <- data_filter %>% filter(logFC < 0)

gsea_all_genes_h <- data$logFC
names(gsea_all_genes_h) <- data$entrezgene_id
gsea_all_genes_h <- gsea_all_genes_h[!is.na(gsea_all_genes_h)]

source("Enrichment.R")

DMD <- enrich_it_baby(gsea_full =  gsea_all_genes_h, 
                       deg_down = data_h_degs_down$entrezgene_id, 
                       deg_up = data_h_degs_up$entrezgene_id, 
                       universe = data$entrezgene_id)



# Étapes 1–2 : Mapping propre
annot <- fData(gset)[, c("ID", "Gene.symbol")]
annot <- annot[!duplicated(annot$ID), ]
expr <- exprs(gset)
expr <- expr[rownames(expr) %in% annot$ID, ]

# Ajouter le symbole de gène comme colonne
expr_df <- as.data.frame(expr)
expr_df$Gene.symbol <- annot$Gene.symbol[match(rownames(expr_df), annot$ID)]

# Étape 3 : Filtrer les symboles valides
expr_df <- expr_df %>%
  filter(!is.na(Gene.symbol), Gene.symbol != "")

# Étape 4 : Agréger par symbole (moyenne)
expr_clean <- expr_df %>%
  group_by(Gene.symbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  column_to_rownames("Gene.symbol") %>%
  as.matrix()

DMD <- list(dataset = data, expr = expr_clean, group = gs)

saveRDS(DMD, file = "papier/data/DMD_paper.RDS")
