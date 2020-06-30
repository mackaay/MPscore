rm(list = ls())

library(limma)
library(RColorBrewer)
pal.grade <- brewer.pal(8, "Dark2")
pal.batch <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Set1")

load("./01TPM_gene.RData")
load("./02eSet_subtype.RData")
load(file = "./sampleinfo_RNA.RData")

sampleinfo$contrast1 <- ifelse(sampleinfo$subtype == "subtype1", "subtype1", "others")
sampleinfo$contrast2 <- ifelse(sampleinfo$subtype == "subtype2", "subtype2", "others")
sampleinfo$contrast3 <- ifelse(sampleinfo$subtype == "subtype3", "subtype3", "others")
sampleinfo$contrast4 <- ifelse(sampleinfo$subtype == "subtype4", "subtype4", "others")



####DEG#####
eSet <-  log2(TPM_gene + 1)
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
DEG1up <- DEG1[which(DEG1$logFC >2 & DEG1$adj.P.Val <0.01),]
DEG1down <- DEG1[which(DEG1$logFC < -2 & DEG1$adj.P.Val <0.01),]

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
DEG2up <- DEG2[which(DEG2$logFC >2 & DEG2$adj.P.Val <0.01),]
DEG2down <- DEG2[which(DEG2$logFC < -2 & DEG2$adj.P.Val <0.01),]


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
DEG3up <- DEG3[which(DEG3$logFC >2 & DEG3$adj.P.Val <0.01),]
DEG3down <- DEG3[which(DEG3$logFC < -2 & DEG3$adj.P.Val <0.01),]


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
DEG4up <- DEG4[which(DEG4$logFC >2 & DEG4$adj.P.Val <0.01),]
DEG4down <- DEG4[which(DEG4$logFC < -2 & DEG4$adj.P.Val <0.01),]

DEG <- c(rownames(DEG1), rownames(DEG2), rownames(DEG3), rownames(DEG4))

DEG <- data.frame(DEG = DEG, subtype = c(rep("subtype1", 357), rep("subtype2",108),
                                         rep("subtype3", 157), rep("subtype4",712)))
DEG$DEG <- as.character(DEG$DEG)
DEG$subtype <- as.character(DEG$subtype)

DEGup <- c(rownames(DEG1up), rownames(DEG2up), rownames(DEG3up), rownames(DEG4up))
DEGup <- data.frame(DEG = DEGup, subtype = c(rep("subtype1", 52), rep("subtype2",80),
                                         rep("subtype3", 16), rep("subtype4",1)), 
                    stringsAsFactors = F)

save(DEGup, file = "./03DEGup_subtype.RData")
save(DEG, file = "./03DEG_subtype.RData")

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
breaksList = seq(-2, 2, by = 0.01)

pdf("./plots/03Heatmap_DEGup.pdf", width = 7, height = 7)
pheatmap(eSet[DEGup$DEG,order(sampleinfo$subtype)], 
         scale = "row",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = F, show_colnames = F,
         fontsize_row = 1, fontsize_col = 1,
         cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         annotation_legend = T,
         annotation_col = annotation_col,
         annotation_colors  = list(
           #Batch = c(GSE85133 =pal.batch[1], GSE136661 =pal.batch[2]),
           Subtype = c(subtype1 = pal.subtype[1], subtype2=pal.subtype[2], 
                       subtype3 = pal.subtype[3], subtype4 = pal.subtype[4]),
           Grade = c(I = pal.grade[1], II =pal.grade[2], III=pal.grade[3])),
         main = "Up-DEG of Meningioma Subtypes"
)
dev.off()
