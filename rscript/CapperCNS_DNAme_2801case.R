rm(list = ls())

library(ff)
library(data.table)

bval <- fread("./GSE90496_Capper/GSE90496_beta.txt.gz", header = T, stringsAsFactors = F)
bval <- as.data.frame(bval)

sampleinfo.MNG <- read.csv("./GSE90496_Capper/sampleinfo_MNG.csv", header = T, stringsAsFactors = F)
colnames(bval)
bval[1:4,1:10]

load("./GSE90496_Capper/sampleinfo_GSE90496_Capper.RData")
sampleinfo <- sampleinfo[grep("MNG" , sampleinfo$`methylation class:ch1`),]

sampleinfo <- sampleinfo[sampleinfo$geo_accession %in% sampleinfo.MNG$id,]

data.frame(sampleinfo$geo_accession, sampleinfo$title, sampleinfo.MNG$id)
sampleinfo.MNG <- sampleinfo.MNG[order(sampleinfo.MNG$id),]
sampleinfo.MNG$sampleID <- sampleinfo$title
sampleinfo.MNG$sampleID <- gsub("MNG, s","S", sampleinfo.MNG$sampleID)
sampleinfo.MNG$sampleID <- substr(sampleinfo.MNG$sampleID, 1L, 11L)
sampleinfo.MNG$sampleID <- gsub(" ","", sampleinfo.MNG$sampleID)
sampleinfo.MNG$sampleID <- gsub("Sample","SAMPLE ", sampleinfo.MNG$sampleID)

bval.MNG <- bval[,sampleinfo.MNG$sampleID ]
rownames(bval.MNG) <- bval$ID_REF

save(bval.MNG, sampleinfo.MNG, file = "./GSE90496_Capper/bval_MNG.RData")
