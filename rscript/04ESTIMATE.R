rm(list = ls())

load("./01TPM_gene.RData")
load("./02eSet_subtype.RData")
load(file ="./Symbol_conversion.RData")
sampleinfo$subtype <- paste("subtype", sampleinfo$subtype, sep = "")

eSet <- as.data.frame(eSet)
eSet <- cbind(rownames(eSet), eSet)
colnames(eSet)[1] <- "symbol"
eSet$symbol <- as.character(eSet$symbol)

write.table(eSet, file = "./ESTIMATE/eSet.txt", sep = "\t", row.names = F)


G_list <- G_list[which(G_list$hgnc_symbol %in% rownames(eSet)),]
G_list <- G_list[!duplicated(G_list$hgnc_symbol),]
eSet <- eSet[rownames(eSet) %in% G_list$hgnc_symbol,]
G_list <- G_list[match(G_list$hgnc_symbol,rownames(eSet)),]

eSet<- eSet[order(rownames(eSet)),]
G_list <- G_list[order(G_list$hgnc_symbol),]
G_list$entrezgene_id <- as.character(G_list$entrezgene_id)
rownames(eSet) <- G_list$entrezgene_id
eSet <- as.data.frame(eSet)
eSet <- cbind(G_list$entrezgene_id, eSet)
colnames(eSet)[1] <- "EntrezID" 
eSet$EntrezID <- as.character(eSet$EntrezID)

write.table(eSet, file = "./ESTIMATE/eSet_entrez.txt", sep = "\t", row.names = F)
test <- eSet[,1:4]
colnames(test) <- paste("test", 1:4, sep = "")
write.table(test, file = "./ESTIMATE/test.txt", sep = "\t", row.names = T)




eSet.file = "./ESTIMATE/eSet.txt"
read.table(eSet.file, sep = "\t", header = T)[1:4,1:4]

keep <- intersect(common_genes$GeneSymbol,rownames(eSet))
library(estimate)
filterCommonGenes(input.f=eSet.file, 
                  output.f="./ESTIMATE/eSet_genes.gct", 
                  id="GeneSymbol")

input.df <- read.table(input.f,
                       header=TRUE,
                       row.names=1,
                       sep="\t", 
                       quote="",
                       stringsAsFactors=FALSE)
rownames(input.df)

estimateScore(input.ds = "./ESTIMATE/eSet_genes.gct",
              output.ds="./ESTIMATE/MNG_estimate_score.gct", 
              platform="illumina")

plotPurity(scores="OV_estimate_score.gct", samples="s516", 
           platform="illumina")

scores=read.table("OV_estimate_score.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
scores
