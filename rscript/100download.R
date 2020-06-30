rm(list = ls())
library(GEOquery)
library(stringr)

GSEset <- getGEO("GSE35069", destdir = "./GSE35069_leukcyte/", getGPL = F)
show(GSEset)
eSet <- exprs(GSEset[[1]])
pData(phenoData(GSEset[[1]]))[,c(2,1)]
colnames(pData(phenoData(GSEset[[1]])))
sampleinfo <- pData(phenoData(GSEset[[1]]))[,c(2,1)]

sampleinfo_GSE35069 <- sampleinfo
bval_GSE35069 <- eSet

save(bval_GSE35069, sampleinfo_GSE35069, file = "./GSE35069_leukcyte/bval.RData")
