library(GEOquery)
GSEset <- getGEO("GSE144037_miR", destdir = "./GSE144037_miR/", getGPL = F)
show(GSEset)
eSet <- exprs(GSEset[[1]])
pData(phenoData(GSEset[[1]]))[,c(2,1)]
sampleinfo <- pData(phenoData(GSEset[[1]]))[,c(2,1)]
sampleinfo$title
sampleinfo$type <- c(rep("retina", 3), rep("RB", 3), rep("invasiveRB",3))
sampleinfo