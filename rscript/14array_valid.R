rm(list = ls())

library(RColorBrewer)
pal.grade <- brewer.pal(8, "Dark2")
pal.batch <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Set1")

load( file = "./03DEG_subtype.RData")
load(file = "./03DEGup_subtype.RData")

DEGup.list <- DEGup$DEG
DEG.list <- DEG$DEG


#####GSE101638 PolyA RNAseq#######




######GSE16181 array ######
#normal meningeal 



######GSE74385 array ######
load("./00GSE74385_array.RData")
eSet.DEGup <- eSet_GSE74385[rownames(eSet_GSE74385) %in% DEGup.list,]
eset.DEG <- eSet_GSE74385[rownames(eSet_GSE74385) %in% DEG.list, ]

library(pheatmap)
col.type <- brewer.pal(9, "Paired")
annotation_col = data.frame(
  Subtype = sampleinfo_GSE74385$`subtype:ch1`
)
rownames(annotation_col) = colnames(eset.DEG)
names(col.type) <- levels(as.factor(sampleinfo_GSE74385$`subtype:ch1`))
breaksList = seq(-2, 2, by = 0.1)

pheatmap(eset.DEG, 
         scale = "row",
         clustering_distance_cols = "manhattan",
         clustering_method = "complete",
         show_rownames = T, show_colnames = F,
         fontsize_row = 10, fontsize_col = 1,
         cluster_cols = T,
         cluster_rows = F,
         annotation_legend = T,
         color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         border_color = NA,
         annotation_col = annotation_col,
         annotation_colors  = list(
           #Time = c("white", "firebrick"),
           # subtype = c(subtype1 = "#A6CEE3", subtype2 = "#1F78B4",
           #            subtype3 = "#B2DF8A", subtype4 = "#33A02C"),
           Subtype = col.type
         ),
         main = "DEGs"
)
