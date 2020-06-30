rm(list = ls())

load("./00GSE84263_array.RData")
load("./03DEG_subtype.RData")

######PCA######
library(limma)
library(RColorBrewer)
pal.celltype <- brewer.pal(8, "Set1")

pdf("./plots/01PCA_GSE84263_array.pdf", width = 5, height = 5)
plotMDS(eSet, top=1000, gene.selection="common", pch = 19, cex = 1.5,
        col=pal.celltype[factor(sampleinfo$tissue)], 
        main = "GSE84263_array")
legend("bottomright", legend=levels(factor(sampleinfo$tissue)), text.col=pal.celltype,
       bg="white", cex=1)
dev.off()

boxplot(eSet)


######subtype cluster######
head(eSet)

library(GEOquery)
gpl <- getGEO('GPL10558', destdir="./GSE84263_array/")
colnames(Table(gpl)) ## [1] 49395 24
head(Table(gpl)[,c(1,13)])
tail(Table(gpl)[,c(1,13)])## you need to check this , which column do you need
probe2symbol=Table(gpl)[,c(1,13)]

probe2symbol <- probe2symbol[which(probe2symbol$ID %in% rownames(eSet)),]
data.frame(probe2symbol$ID, rownames(eSet))
rownames(eSet)<- probe2symbol$Symbol

eSet_array <- eSet[which(rownames(eSet) %in% DEG$DEG),]


library(RColorBrewer)
pal.grade <- brewer.pal(8, "Dark2")
pal.batch <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Set1")

library(pheatmap)

annotation_col = data.frame(
  Grade = sampleinfo$grade,
  Batch = sampleinfo$batch,
  Subtype = sampleinfo$subtype
)
rownames(annotation_col) = colnames(eSet)


pdf("./plots/03Heatmap_subtypeDEG_array.pdf", width = 7, height = 7)
breaksList = seq(-3, 3, by = 0.01)
pheatmap(eSet_array, 
         scale = "row",
         #kmeans_k = 4,
         border_color = NA,
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D",
         show_rownames = F, show_colnames = F,
         fontsize_row = 1, fontsize_col = 1,
         cluster_cols = T,
         cluster_rows = T,
         color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
#         annotation_legend = T,
#         annotation_col = annotation_col,
#         annotation_colors  = list(
#           Batch = c(GSE85133 =pal.batch[1], GSE136661 =pal.batch[2]),
#           Subtype = c(subtype1 = pal.subtype[1], subtype2=pal.subtype[2], 
#                       subtype3 = pal.subtype[3], subtype4 = pal.subtype[4]),
#           Grade = c(I = pal.grade[1], II =pal.grade[2], III=pal.grade[3])),
         main = "Meningioma (microarray)"
)
