rm(list = ls())

load("./GSE90496_Capper/bval_MNG.RData")
library(limma)
library(RColorBrewer)
library(ff)
library(data.table)

pal.grade <- brewer.pal(8, "Dark2")


library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

LUMPcg <- read.csv("./LUMPcg.csv", header = T, stringsAsFactors = F)

keep <- LUMPcg[LUMPcg$LUMP %in% rownames(ann450k),]

######LeMoN score######
LUMP.bval <- bval.MNG[which(rownames(bval.MNG) %in% keep),]
rowMeans(LUMP.bval, na.rm = T)
LeMoN.cg <- rownames(LUMP.bval)

bval.MNG.lemon <- bval.MNG[LeMoN.cg , ]

min(colMeans(bval.MNG.lemon, na.rm = T)/0.85,1)



####heatmap of LeMoN score cg 
load("./GSE35069_leukcyte/bval.RData")

bval_GSE35069 <- bval_GSE35069[LeMoN.cg,]
rownames(bval_GSE35069)
rownames(bval.MNG.lemon)

bval.heatmap <- cbind(bval.MNG.lemon, bval_GSE35069)
colnames(bval.heatmap)
sample.ann <- data.frame(ID = colnames(bval.heatmap), 
                         type = c(rep("MNG", 90), rep("Leukocyte", 60)))
rownames(sample.ann) <- sample.ann$ID


pal.type <- brewer.pal(8, "Set1")

library(pheatmap)

annotation_col = data.frame(
  type = sample.ann$type
)
rownames(annotation_col) = colnames(bval.heatmap)

pdf("./plots/102Heatmap_LeMoN.pdf", width = 5, height = 5)
breaksList = seq(-3, 3, by = 0.01)
pheatmap(bval.heatmap, 
         #scale = "row",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = F, show_colnames = F,
         fontsize_row = 1, fontsize_col = 1,
         cluster_cols = F,
         cluster_rows = T,
         border_color = NA,
         color = colorRampPalette(c("darkblue", "black", "yellow"))(50), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         #breaks = breaksList,
         annotation_legend = T,
         annotation_col = annotation_col,
         annotation_colors  = list(
           type = c(MNG =pal.type[1], Leukocyte =pal.type[2])
           ),
         main = ""
)
dev.off()



#####correlation LeMoN with tumor purity