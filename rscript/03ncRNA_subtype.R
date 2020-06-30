rm(list = ls())

library(limma)

load("./01TPM_gene.RData")

load("./02eSet_subtype.RData")

sampleinfo$subtype <- paste("subtype", sampleinfo$subtype, sep = "")

######lncRNA#####
eSet <- log2(TPM_lncRNA+1)

sampleinfo$contrast1 <- ifelse(sampleinfo$subtype == "subtype1", "subtype1", "others")
sampleinfo$contrast2 <- ifelse(sampleinfo$subtype == "subtype2", "subtype2", "others")
sampleinfo$contrast3 <- ifelse(sampleinfo$subtype == "subtype3", "subtype3", "others")
sampleinfo$contrast4 <- ifelse(sampleinfo$subtype == "subtype4", "subtype4", "others")


cellType <- factor(sampleinfo$contrast1)
table(cellType)

design <- model.matrix(~0+cellType, data=sampleinfo)
colnames(design) <- c(levels(cellType))

fit <- lmFit(eSet, design)

contMatrix <- makeContrasts(subtype1-others,
                            levels=design)
contMatrix

fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))

DEG1 <- topTable(fit2, num=Inf, coef=1)
DEG1 <- DEG1[which(abs(DEG1$logFC)>1.5 & DEG1$adj.P.Val <0.05),]


cellType <- factor(sampleinfo$contrast2)
table(cellType)

design <- model.matrix(~0+cellType, data=sampleinfo)
colnames(design) <- c(levels(cellType))

fit <- lmFit(eSet, design)

contMatrix <- makeContrasts(subtype2-others,
                            levels=design)
contMatrix

fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))

DEG2 <- topTable(fit2, num=Inf, coef=1)
DEG2 <- DEG2[which(abs(DEG2$logFC)>1.5 & DEG2$adj.P.Val <0.05),]



cellType <- factor(sampleinfo$contrast3)
table(cellType)

design <- model.matrix(~0+cellType, data=sampleinfo)
colnames(design) <- c(levels(cellType))

fit <- lmFit(eSet, design)

contMatrix <- makeContrasts(subtype3-others,
                            levels=design)
contMatrix

fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))

DEG3 <- topTable(fit2, num=Inf, coef=1)
DEG3 <- DEG3[which(abs(DEG3$logFC)>1.5 & DEG3$adj.P.Val <0.05),]



cellType <- factor(sampleinfo$contrast4)
table(cellType)

design <- model.matrix(~0+cellType, data=sampleinfo)
colnames(design) <- c(levels(cellType))

fit <- lmFit(eSet, design)

contMatrix <- makeContrasts(subtype4-others,
                            levels=design)
contMatrix

fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))

DEG4 <- topTable(fit2, num=Inf, coef=1)
DEG4 <- DEG4[which(abs(DEG4$logFC)>1.5 & DEG4$adj.P.Val <0.05),]

#DEG1<- DEG1[which(DEG1$logFC >0),]
#DEG2<- DEG2[which(DEG2$logFC >0),]
#DEG3<- DEG3[which(DEG3$logFC >0),]
#DEG4<- DEG4[which(DEG4$logFC >0),]
DEG <- c(rownames(DEG1), rownames(DEG2), rownames(DEG3), rownames(DEG4))

DEG <- data.frame(DEG = DEG, subtype = c(rep("subtype1", nrow(DEG1)), rep("subtype2",nrow(DEG2)),
                                         rep("subtype3", nrow(DEG3)), rep("subtype4",nrow(DEG4))))
DEG$DEG <- as.character(DEG$DEG)
DEG$subtype <- as.character(DEG$subtype)

DElncRNA <- DEG

save(DElncRNA, file = "03DElncRNA_subtype.RData")
write.csv(DElncRNA, file = "./DElncRNA.csv")
#######heatmap ######
library(RColorBrewer)
pal.grade <- brewer.pal(8, "Dark2")
pal.batch <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Set1")


library(pheatmap)

annotation_col = data.frame(
  Grade = sampleinfo$grade,
  #Batch = sampleinfo$batch,
  Subtype = sampleinfo$subtype
)
rownames(annotation_col) = colnames(eSet)

pdf("./plots/03Heatmap_subtypelncRNA.pdf", width = 5, height = 5)
breaksList = seq(-4, 4, by = 0.01)
pheatmap(eSet[DEG$DEG,order(sampleinfo$subtype)], 
         scale = "row",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = F, show_colnames = F,
         fontsize_row = 1, fontsize_col = 1,
         cluster_cols = F,
         cluster_rows = T,
         color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         annotation_legend = T,
         annotation_col = annotation_col,
         annotation_colors  = list(
           #Batch = c(GSE85133 =pal.batch[1], GSE136661 =pal.batch[2]),
           Subtype = c(subtype1 = pal.subtype[1], subtype2=pal.subtype[2], 
                       subtype3 = pal.subtype[3], subtype4 = pal.subtype[4]),
           Grade = c(I = pal.grade[1], II =pal.grade[2], III=pal.grade[3])),
         main = "Meningioma Subtypes lncRNA"
)
dev.off()


######miRNA#####
eSet <- log2(TPM_miRNA+1)

sampleinfo$contrast1 <- ifelse(sampleinfo$subtype == "subtype1", "subtype1", "others")
sampleinfo$contrast2 <- ifelse(sampleinfo$subtype == "subtype2", "subtype2", "others")
sampleinfo$contrast3 <- ifelse(sampleinfo$subtype == "subtype3", "subtype3", "others")
sampleinfo$contrast4 <- ifelse(sampleinfo$subtype == "subtype4", "subtype4", "others")


cellType <- factor(sampleinfo$contrast1)
table(cellType)

design <- model.matrix(~0+cellType, data=sampleinfo)
colnames(design) <- c(levels(cellType))

fit <- lmFit(eSet, design)

contMatrix <- makeContrasts(subtype1-others,
                            levels=design)
contMatrix

fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))

DEG1 <- topTable(fit2, num=Inf, coef=1)
DEG1 <- DEG1[which(abs(DEG1$logFC)>1 & DEG1$adj.P.Val <0.05),]


cellType <- factor(sampleinfo$contrast2)
table(cellType)

design <- model.matrix(~0+cellType, data=sampleinfo)
colnames(design) <- c(levels(cellType))

fit <- lmFit(eSet, design)

contMatrix <- makeContrasts(subtype2-others,
                            levels=design)
contMatrix

fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))

DEG2 <- topTable(fit2, num=Inf, coef=1)
DEG2 <- DEG2[which(abs(DEG2$logFC)>1 & DEG2$adj.P.Val <0.05),]



cellType <- factor(sampleinfo$contrast3)
table(cellType)

design <- model.matrix(~0+cellType, data=sampleinfo)
colnames(design) <- c(levels(cellType))

fit <- lmFit(eSet, design)

contMatrix <- makeContrasts(subtype3-others,
                            levels=design)
contMatrix

fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))

DEG3 <- topTable(fit2, num=Inf, coef=1)
DEG3 <- DEG3[which(abs(DEG3$logFC)>1 & DEG3$adj.P.Val <0.05),]



cellType <- factor(sampleinfo$contrast4)
table(cellType)

design <- model.matrix(~0+cellType, data=sampleinfo)
colnames(design) <- c(levels(cellType))

fit <- lmFit(eSet, design)

contMatrix <- makeContrasts(subtype4-others,
                            levels=design)
contMatrix

fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))

DEG4 <- topTable(fit2, num=Inf, coef=1)
DEG4 <- DEG4[which(abs(DEG4$logFC)>1 & DEG4$adj.P.Val <0.05),]

#DEG1<- DEG1[which(DEG1$logFC >0),]
#DEG2<- DEG2[which(DEG2$logFC >0),]
#DEG3<- DEG3[which(DEG3$logFC >0),]
#DEG4<- DEG4[which(DEG4$logFC >0),]
DEG <- c(rownames(DEG1), rownames(DEG2), rownames(DEG3), rownames(DEG4))

DEG <- data.frame(DEG = DEG, subtype = c(rep("subtype1", nrow(DEG1)), rep("subtype2",nrow(DEG2)),
                                         rep("subtype3", nrow(DEG3)), rep("subtype4",nrow(DEG4))))
DEG$DEG <- as.character(DEG$DEG)
DEG$subtype <- as.character(DEG$subtype)

DEmiRNA <- DEG

save(DEmiRNA, file = "03DEmiRNA_subtype.RData")

#######heatmap ######
library(RColorBrewer)
pal.grade <- brewer.pal(8, "Dark2")
pal.batch <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Set1")


library(pheatmap)

annotation_col = data.frame(
  Grade = sampleinfo$grade,
  #Batch = sampleinfo$batch,
  Subtype = sampleinfo$subtype
)
rownames(annotation_col) = colnames(eSet)

pdf("./plots/03Heatmap_subtypemiRNA.pdf", width = 5, height = 5)
breaksList = seq(-4, 4, by = 0.01)
pheatmap(eSet[DEG$DEG,order(sampleinfo$subtype)], 
         scale = "row",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = F, show_colnames = F,
         fontsize_row = 1, fontsize_col = 1,
         cluster_cols = F,
         cluster_rows = T,
         color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         annotation_legend = T,
         annotation_col = annotation_col,
         annotation_colors  = list(
           #Batch = c(GSE85133 =pal.batch[1], GSE136661 =pal.batch[2]),
           Subtype = c(subtype1 = pal.subtype[1], subtype2=pal.subtype[2], 
                       subtype3 = pal.subtype[3], subtype4 = pal.subtype[4]),
           Grade = c(I = pal.grade[1], II =pal.grade[2], III=pal.grade[3])),
         main = "Meningioma Subtypes miRNA"
)
dev.off()



######unsupervised variable clustering######
library(pheatmap)
library(RColorBrewer)
pal.grade <- brewer.pal(8, "Dark2")
pal.batch <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Set1")
##   lncRNA 
nrow(TPM_lncRNA)
nrow(TPM_miRNA)

eSet <- log2(TPM_lncRNA+1)
eSet <- eSet[rowMeans(eSet)> 0,] #mean of TPM larger than 1
highly.variable <- apply(eSet, 1, sd, na.rm = TRUE)
highly.variable <- sort(highly.variable, decreasing = T)
highly.variable <- highly.variable[1:round(0.01*length(highly.variable))]
eSet <- eSet[names(highly.variable),]

annotation_col = data.frame(
  Grade = sampleinfo$grade,
  #Batch = sampleinfo$batch,
  Subtype = sampleinfo$subtype
)
rownames(annotation_col) = colnames(eSet)

breaksList = seq(-4, 4, by = 0.01)
pdf("./plots/03Heatmap_lncRNA_unsupervised.pdf", width = 5, height = 5)
pheatmap(eSet, 
         scale = "row",
         #kmeans_k = 4,
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D",
         show_rownames = T, show_colnames = F,
         fontsize_row = 5, fontsize_col = 1,
         cluster_cols = T,
         cluster_rows = T,
         color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         annotation_legend = T,
         annotation_col = annotation_col,
         annotation_colors  = list(
           #Batch = c(GSE85133 =pal.batch[1], GSE136661 =pal.batch[2]),
           Subtype = c(subtype1 = pal.subtype[1], subtype2=pal.subtype[2], 
                       subtype3 = pal.subtype[3], subtype4 = pal.subtype[4]),
           Grade = c(I = pal.grade[1], II =pal.grade[2], III=pal.grade[3])),
         main = "Meningioma Subtypes lncRNA"
)
dev.off()



### miRNA
eSet <- log2(TPM_miRNA+1)
eSet <- eSet[rowMeans(eSet)> 0,] #mean of TPM larger than 1
highly.variable <- apply(eSet, 1, sd, na.rm = TRUE)
highly.variable <- sort(highly.variable, decreasing = T)
highly.variable <- highly.variable[1:round(0.1*length(highly.variable))]
eSet <- eSet[names(highly.variable),]

annotation_col = data.frame(
  Grade = sampleinfo$grade,
  #Batch = sampleinfo$batch,
  Subtype = sampleinfo$subtype
)
rownames(annotation_col) = colnames(eSet)

breaksList = seq(-4, 4, by = 0.01)
pdf("./plots/03Heatmap_miRNA_unsupervised.pdf", width = 5, height = 5)
pheatmap(eSet, 
         scale = "row",
         #kmeans_k = 4,
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D",
         show_rownames = T, show_colnames = F,
         fontsize_row = 5, fontsize_col = 1,
         cluster_cols = T,
         cluster_rows = T,
         color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         annotation_legend = T,
         annotation_col = annotation_col,
         annotation_colors  = list(
           #Batch = c(GSE85133 =pal.batch[1], GSE136661 =pal.batch[2]),
           Subtype = c(subtype1 = pal.subtype[1], subtype2=pal.subtype[2], 
                       subtype3 = pal.subtype[3], subtype4 = pal.subtype[4]),
           Grade = c(I = pal.grade[1], II =pal.grade[2], III=pal.grade[3])),
         main = "Meningioma Subtypes miRNA"
)
dev.off()


