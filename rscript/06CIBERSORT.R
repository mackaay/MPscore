rm(list = ls())

#library(limma)
library(RColorBrewer)
pal.grade <- brewer.pal(8, "Dark2")
pal.batch <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Set1")

load(file = "./TPM.RData")
load("./01TPM_gene.RData")
load("./02eSet_subtype.RData")
load(file = "./sampleinfo_RNA.RData")

tmp <- read.csv("./CHOL_TCGA.csv", sep = ",", header = T, stringsAsFactors = F)
symbol <- tmp$SYMBOL

TPM_gene <- TPM[rownames(TPM) %in% symbol,]

write.csv(TPM_gene, file = "./TPM_gene.csv", row.names = T)
save(TPM_gene, file = "./06TPM_gene.RData")

setdiff(tmp$SYMBOL, rownames(TPM_gene))
grep("BCORP1", rownames(TPM))




#####matrix download from GEO######
load(file = "./00GSE136661_RNAseq.RData")

rownames(counts)
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensg <- rownames(counts)
symbol <- c()
G_list <- getBM(filters= "ensembl_gene_id", 
                attributes= c("ensembl_gene_id", 
                              "hgnc_symbol", 
                              "entrezgene_id"),
                values=ensg,
                mart= mart)
save(G_list, file ="./Symbol_conversion.RData")

counts <- counts[rownames(counts) %in% G_list$ensembl_gene_id,]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
G_list <- G_list[match(G_list$ensembl_gene_id, rownames(counts)),]

setdiff(rownames(counts), G_list$ensembl_gene_id)
counts <- counts[rownames(counts) %in% G_list$ensembl_gene_id,]
counts$symbol <- G_list$hgnc_symbol
counts <- counts[!duplicated(counts$symbol),]
rownames(counts) <- counts$symbol
counts <- counts[,1:160]

save(counts, sampleinfo, file = "./00GSE136661_RNAseq_GEOmatrix.RData")

tmp <- read.table("./CHOL_TCGA.txt", sep = "\t", header = T, stringsAsFactors = F)
symbol <- tmp$SYMBOL

counts <- counts[rownames(counts) %in% symbol,]

write.table(counts, file = "./GSE136661_GEOmatrix.txt", row.names = T, sep = "\t")
