rm(list = ls())


rawdata <- read.delim("./GSE136661_RNAseq_raw/RNA_name.txt", header = F, sep = "\t",
                      stringsAsFactors = F)
tmp <- read.delim("./GSE85133_RNAseq_raw/RNA_name.txt", header = F, sep = "\t",
                      stringsAsFactors = F)

rawdata <- cbind(tmp, rawdata[7:166])
rawdata <- rawdata[-1,]

colnames(rawdata) <- rawdata[1,]

rownames(rawdata) <- rawdata$Geneid

colnames(rawdata)[7:185] <- substr(colnames(rawdata)[7:185], 13,23)
colnames(rawdata)
colnames(rawdata)[7:185] <- gsub("_", "", colnames(rawdata)[7:185])

rawdata <- rawdata[-1,]
gene_kilobase <- as.numeric(rawdata$Length)/1000
tmp <- as.data.frame(rawdata[,7:185])
tmp <- apply(tmp, 2, as.numeric)
rownames(tmp) <- rawdata$Geneid

RPK <- tmp/gene_kilobase
coverage <- colSums(RPK[,1:179])/1000000
RPK <- t(RPK)
TPM <- t(RPK/coverage)
TPM <- round(TPM, 4)
colSums(TPM)
head(TPM)

save(TPM, file = "./TPM.RData")
