rm(list = ls())
library(GEOquery)
library(stringr)

######GSE85133 RNA seq######
GSEset <- getGEO("GSE85133", destdir = "./GSE85133_RNAseq/", getGPL = F)
show(GSEset)
eSet <- exprs(GSEset[[1]])
pData(phenoData(GSEset[[1]]))[,c(2,1)]
colnames(pData(phenoData(GSEset[[1]])))
sampleinfo <- pData(phenoData(GSEset[[1]]))[,c(2,1, 44:50)]
sampleinfo$title
sampleinfo


dir <- './GSE85133_RNAseq/'

FPKM <- read.table("./GSE85133_RNAseq/GSM2258207_MN-1037_TUMOR_genes.fpkm_tracking.txt.gz", 
                   sep = '\t', header = T)
gene_id <-  FPKM$gene_id

FPKM <- as.data.frame(FPKM$FPKM)
colnames(FPKM) <-  str_sub(basename(".AnaM/GSE85133_RNAseq/GSM2258207_MN-1037_TUMOR_genes.fpkm_tracking.txt.gz")
                           , start = 1L, end = 10L)
FPKM$gene_id <- gene_id

head(FPKM)

for (i in list.files(path = dir,pattern = '*ing.txt.gz$', recursive = T)[2:19]) {
  tmp2 <- read.table(paste(dir,i,sep = ""), 
                     sep = '\t', header = T)
  tmp2 <- as.data.frame(tmp2$FPKM)
  colnames(tmp2) <-  str_sub(basename(i)
                             , start = 1L, end = 10L)
  FPKM <- cbind(FPKM, tmp2)
  print(i)
}

sampleinfo_GSE85133 <- sampleinfo
save(sampleinfo_GSE85133, file = "./00GSE85133_RNAseq.RData" )


######GSE136661 RNA seq######
GSEset <- getGEO("GSE136661", destdir = "./GSE136661_RNAseq/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
colnames(pData(phenoData(GSEset[[1]])))
sampleinfo <- pData(phenoData(GSEset[[1]]))[,c(2,1, 51:56)]
sampleinfo
colnames(sampleinfo)[3:8] <- c("age", "cohort", "gender", "pathology", 
                               "tissue", "wes")


dir <- './GSE136661_RNAseq/'

counts <- read.table("./GSE136661_RNAseq/GSM4054837_12-25748.20R_htseq.counts.txt.gz", 
                   sep = '\t', header = F)

ensembl <-  counts$V1

counts <- as.data.frame(counts$V2)
colnames(counts) <-  str_sub(basename("./GSE136661_RNAseq/GSM4054837_12-25748.20R_htseq.counts.txt.gz"), 
                           start = 1L, end = 10L)
rownames(counts) <- ensembl
head(counts)

for (i in list.files(path = dir,pattern = '*counts.txt.gz$', recursive = T)[2:160]) {
  tmp2 <- read.table(paste(dir,i,sep = ""), 
                     sep = '\t', header = F)
  tmp2 <- as.data.frame(tmp2$V2)
  colnames(tmp2) <-  str_sub(basename(i),
                             start = 1L, end = 10L)
  counts <- cbind(counts, tmp2)
  print(i)
}

save(counts, sampleinfo, file = "./00GSE136661_RNAseq.RData")



######GSE84263 array ######
GSEset <- getGEO("GSE84263", destdir = "./GSE84263_array/", getGPL = F)
show(GSEset)
eSet <- exprs(GSEset[[1]])
pData(phenoData(GSEset[[1]]))[,c(2,1)]
colnames(pData(phenoData(GSEset[[1]])))
sampleinfo <- pData(phenoData(GSEset[[1]]))[,c(2,1, 38:43)]
sampleinfo
colnames(sampleinfo)[3:8] <- c("batch", "grade", "mutation", 
                               "primary", "radiation", "tissue")

library("biomaRt")
ensembl <- useMart("ensembl")
head(listDatasets(ensembl))
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
View(listAttributes(ensembl))
probe2symbol <- getBM(attributes = c("illumina_humanht_12_v4", "hgnc_symbol"), 
                      filters = "illumina_humanht_12_v4", 
                      values = rownames(eSet), 
                      mart = ensembl)
probe.ov <- intersect(rownames(eSet) , probe2symbol$illumina_humanht_12_v4)
eSet <- eSet[probe.ov,]
probe2symbol <- probe2symbol[probe2symbol$illumina_humanht_12_v4 %in% probe.ov,]
probe2symbol <- probe2symbol[match(rownames(eSet), probe2symbol$illumina_humanht_12_v4),]
eSet <- as.data.frame(eSet)
eSet$symbol <- probe2symbol$hgnc_symbol
eSet <- eSet[!duplicated(eSet$symbol),]
eSet <- eSet[,-ncol(eSet)]

save(eSet, sampleinfo, file = "./00GSE84263_array.RData")




######GSE16581 array ######
GSEset <- getGEO("GSE16581", destdir = "./GSE16581_array/", getGPL = F)
show(GSEset)
eSet <- exprs(GSEset[[1]])
pData(phenoData(GSEset[[1]]))[,c(2,1:16)]
colnames(pData(phenoData(GSEset[[1]])))
sampleinfo <- pData(phenoData(GSEset[[1]]))[,c(2,1, 44:50)]
sampleinfo

gpl <- getGEO('GPL570', destdir="./GSE16581_array/")
colnames(Table(gpl)) ## [1] 49395 24
head(Table(gpl)[,c(1,11)]) ## you need to check this , which column do you need
head(Table(gpl)[1:4,])
probe2symbol=Table(gpl)[,c(1,11)]

eSet <- eSet[match(rownames(eSet),probe2symbol$ID),]
eSet <- as.data.frame(eSet)
eSet$symbol <- probe2symbol$`Gene Symbol`

save(eSet, sampleinfo, file = "./00GSE16581_array.RData")


######GSE16181 array ######
GSEset <- getGEO("GSE16181", destdir = "./GSE16181_array/", getGPL = F)
show(GSEset)
eSet <- exprs(GSEset[[1]])
pData(phenoData(GSEset[[1]]))[,c(2,1:16)]
colnames(pData(phenoData(GSEset[[1]])))
sampleinfo <- pData(phenoData(GSEset[[1]]))[,c(2,1, 47:55)]
sampleinfo

gpl <- getGEO('GPL8557', destdir="./GSE16181_array/")
colnames(Table(gpl)) ## [1] 49395 24
head(Table(gpl)[,c(1,2)]) ## you need to check this , which column do you need
head(Table(gpl)[1:4,])
probe2symbol=Table(gpl)[,c(1,2)]

eSet <- eSet[match(rownames(eSet),probe2symbol$ID),]
eSet <- as.data.frame(eSet)
rownames(eSet) <- probe2symbol$GeneName

eSet_GSE16181 <- eSet
sampleinfo_GSE16181 <- sampleinfo

save(eSet_GSE16181, sampleinfo_GSE16181, file = "./00GSE16181_array.RData")


######GSE43290 array ######
# this info from GEO is missing
GSEset <- getGEO("GSE43290", destdir = "./GSE43290_array/", getGPL = F)
show(GSEset)
eSet <- exprs(GSEset[[1]])
pData(phenoData(GSEset[[1]]))[,c(2,1:16)]
colnames(pData(phenoData(GSEset[[1]])))
sampleinfo <- pData(phenoData(GSEset[[1]]))[,c(2,1, 42:49)]
sampleinfo

gpl <- getGEO('GPL96', destdir="./GSE43290_array/")
colnames(Table(gpl)) ## [1] 49395 24
head(Table(gpl)[,c(1,2)]) ## you need to check this , which column do you need
head(Table(gpl)[1:4,])
probe2symbol=Table(gpl)[,c(1,2)]

eSet <- eSet[match(rownames(eSet),probe2symbol$ID),]
eSet <- as.data.frame(eSet)
rownames(eSet) <- probe2symbol$GeneName

eSet_GSE16181 <- eSet
sampleinfo_GSE16181 <- sampleinfo

save(eSet_GSE16181, sampleinfo_GSE16181, file = "./00GSE43290_array.RData")


######GSE74385_array ######
GSEset <- getGEO("GSE74385", destdir = "./GSE74385_array/", getGPL = F)
show(GSEset)
eSet <- exprs(GSEset[[1]])
pData(phenoData(GSEset[[1]]))[,c(2,1:16)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
sampleinfo <- pData(phenoData(GSEset[[1]]))[,c(2,1, 33)]
sampleinfo

library("biomaRt")
ensembl <- useMart("ensembl")
head(listDatasets(ensembl))
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
View(listAttributes(ensembl))
probe2symbol <- getBM(attributes = c("illumina_humanht_12_v4", "hgnc_symbol"), 
                      filters = "illumina_humanht_12_v4", 
                      values = rownames(eSet), 
                      mart = ensembl)
probe.ov <- intersect(rownames(eSet) , probe2symbol$illumina_humanht_12_v4)
eSet <- eSet[probe.ov,]
eSet <- as.data.frame(eSet)
probe2symbol <- probe2symbol[!duplicated(probe2symbol$illumina_humanht_12_v4),]
rownames(probe2symbol) <- probe2symbol$illumina_humanht_12_v4
probe2symbol <- probe2symbol[ probe.ov ,]
eSet <- eSet[match(rownames(eSet),probe2symbol$illumina_humanht_12_v4),]

eSet$symbol <- probe2symbol$hgnc_symbol
eSet <- eSet[!duplicated(eSet$symbol),]
rownames(eSet) <- eSet$symbol
eSet  <- eSet[,-ncol(eSet)]

eSet_GSE74385 <- eSet
sampleinfo_GSE74385 <- sampleinfo

save(eSet_GSE74385, sampleinfo_GSE74385, file = "./00GSE74385_array.RData")


######GSE101638 PolyA RNA seq ######
GSEset <- getGEO("GSE101638", destdir = "./GSE101638_PolyA_RNAseq/", getGPL = F)
show(GSEset)
#eSet <- exprs(GSEset[[1]])
pData(phenoData(GSEset[[1]]))[,c(2,1:16)]
colnames(pData(phenoData(GSEset[[1]])))
sampleinfo <- pData(phenoData(GSEset[[1]]))[,c(2,1)]
sampleinfo

write.csv(sampleinfo, file = "./GSE101638_PolyA_RNAseq/sampleinfo.csv", row.names = T)


######GSE74385 ######
GSEset <- getGEO("GSE74385", destdir = "./GSE74385_validation/", getGPL = F)
show(GSEset)
eSet <- exprs(GSEset[[1]])
pData(phenoData(GSEset[[1]]))[,c(2,1,33)]
colnames(pData(phenoData(GSEset[[1]])))
sampleinfo <- pData(phenoData(GSEset[[1]]))[,c(2,1, 33)]
sampleinfo

gpl <- getGEO('GPL10558', destdir="./GSE74385_validation/")
colnames(Table(gpl)) ## [1] 49395 24
head(Table(gpl)[,c(1,13)]) ## you need to check this , which column do you need
head(Table(gpl)[1:4,])
probe2symbol=Table(gpl)[,c(1,13)]

probe2symbol <- probe2symbol[probe2symbol$ID %in% rownames(eSet),]
eSet <- eSet[match(rownames(eSet), probe2symbol$ID),]
rownames(eSet) <- probe2symbol$Symbol

save(eSet, sampleinfo, file = "./00GSE74385_validation.RData")
