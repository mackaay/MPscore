rm(list = ls())

library(limma)
library(RColorBrewer)

load("./TPM.RData")
load(file = "./00GSE136661_RNAseq.RData")
load(file = "./00GSE85133_RNAseq.RData")



#match sampleinfo and TPM
sampleinfo <- data.frame(accession = sampleinfo$geo_accession, 
                         age = sampleinfo$age, 
                         title = sampleinfo$title, 
                         gender = sampleinfo$gender,
                         tissue = sampleinfo$tissue, 
                         grade = sampleinfo$pathology,
                         mutation = rep(NA, 160))
sampleinfo_GSE85133 <- data.frame(accession = sampleinfo_GSE85133$geo_accession, 
                         age = sampleinfo_GSE85133$`age:ch1`, 
                         title = sampleinfo_GSE85133$title, 
                         gender = sampleinfo_GSE85133$`gender:ch1`,
                         tissue = sampleinfo_GSE85133$`tissue:ch1`, 
                         grade = sampleinfo_GSE85133$`who grade:ch1`,
                         mutation = sampleinfo_GSE85133$`mutation group:ch1`)
sampleinfo <- rbind(sampleinfo_GSE85133, sampleinfo)


sratable <- read.delim("./GSE85133_RNAseq_raw/SraRunTable.txt", header = T, sep = ",",
                      stringsAsFactors = F)
sratable <- cbind(sratable$Run, sratable$GEO_Accession..exp.)
tmp <- read.delim("./GSE136661_RNAseq_raw/SraRunTable.txt", header = T, sep = ",",
                       stringsAsFactors = F)
tmp <- cbind(tmp$Run, tmp$GEO_Accession..exp.)

sratable <- rbind(sratable, tmp)
colnames(sratable) <- c("sra", "geo")

sratable <- sratable[order(as.numeric(gsub("SRR","",sratable[,1]))),]

data.frame(colnames(TPM), sratable[,1], sratable[,2] , sampleinfo$accession)


sampleinfo$sra <- colnames(TPM)
rownames(sampleinfo) <- sampleinfo$sra


sampleinfo$gender <- gsub("Female", "F", sampleinfo$gender)
sampleinfo$gender <- gsub("Male", "M", sampleinfo$gender)
sampleinfo$grade <- gsub("1", "I", sampleinfo$grade)
sampleinfo$grade <- gsub("WHO I", "I", sampleinfo$grade)
sampleinfo$batch <- c(rep("GSE85133",19), rep("GSE136661",160))


#gene list
gene_namelist <- read.delim("./protein-coding_gene.txt", 
                            header = T, sep = "\t", stringsAsFactors = F)
lncRNA_list <- read.delim("./RNA_long_non-coding.txt", 
                          header = T, sep = "\t", stringsAsFactors = F)
miRNA_list <- read.delim("./RNA_micro.txt", 
                         header = T, sep = "\t", stringsAsFactors = F)

TPM_gene <- TPM[which(rownames(TPM) %in% gene_namelist$symbol),]
TPM_lncRNA <- TPM[which(rownames(TPM) %in% lncRNA_list$symbol),]
TPM_miRNA <- TPM[which(rownames(TPM) %in% miRNA_list$symbol),]


#normalization 
TPM_gene.rmBatch <- removeBatchEffect(log2(TPM_gene+0.0001), sampleinfo$batch)

#check 
pdf("./plots/01Boxplot_raw.pdf", width = 10, height = 5)
par(cex.lab=1.5) # is for y-axis
par(cex.axis=0.5) # is for x-axis
boxplot(log2(TPM_gene+0.0001), main = "Meningioma Protein Coding Genes", ylab = "log2(TPM)",las=2)
boxplot(TPM_gene.rmBatch, main = "Meningioma Protein Coding Genes", ylab = "log2(TPM)",las=2)
boxplot(log2(TPM_lncRNA+0.0001), main = "Meningioma lncRNAs", ylab = "log2(TPM)",las=2)
boxplot(log2(TPM_miRNA+0.0001), main = "Meningioma miRNAs", ylab = "log2(TPM)",las=2)
dev.off()

pal.grade <- brewer.pal(8, "Dark2")
pal.batch <- brewer.pal(8, "Accent")
pal.gender <- brewer.pal(8, "Set1")

pdf("./plots/01PCA.pdf", width = 5, height = 5)
plotMDS(TPM_gene.rmBatch, top=1000, gene.selection="common", 
        col=pal.grade[factor(sampleinfo$grade)], 
        pch = 19, cex = 1.2,
        main = "WHO Grade")
legend("topleft", legend=levels(factor(sampleinfo$grade)), text.col=pal.grade,
       bg="white", cex=1)

plotMDS(TPM_gene.rmBatch, top=1000, gene.selection="common", 
        col=pal.batch[factor(sampleinfo$batch)], 
        pch = 19, cex = 1.2,
        main = "Batch")
legend("topleft", legend=levels(factor(sampleinfo$batch)), text.col=pal.batch,
       bg="white", cex=1)
dev.off()


save(pal.grade, pal.batch, file = "./color.RData")
save(TPM_gene, TPM_lncRNA, TPM_miRNA, sampleinfo, file = "./01TPM_gene.RData")


