rm(list = ls())

library(limma)
library(RColorBrewer)
library(ff)
library(data.table)

pal.grade <- brewer.pal(8, "Dark2")
pal.type <- brewer.pal(8, "Set1")
pal.grade <- brewer.pal(8, "Accent")

#### CNS methylation
load("./GSE90496_Capper/sampleinfo_GSE90496_Capper.RData")
which(!sampleinfo$`methylation class:ch1` %in% grep("CONTR", sampleinfo$`methylation class:ch1`))

bval <- fread("./GSE90496_Capper/GSE90496_beta.txt.gz", header = T, stringsAsFactors = F)
bval <- as.data.frame(bval)
rownames(bval) <- bval$ID_REF
bval <- bval[,grep("SAMPLE", colnames(bval))]

#bval <- bval[,which(!sampleinfo$`methylation class:ch1` %in% grep("CONTR", sampleinfo$`methylation class:ch1`))]


CNS.me <- data.frame(bval.min = apply(bval, 1, mean, na.rm=T),
                     CNS.me = NA,
                     cg = rownames(bval), 
                     stringsAsFactors = F)
CNS.me$CNS.me <- ifelse(CNS.me$bval.min > 0.5, "LUMP", "No")
CNS.me <- CNS.me[CNS.me$CNS.me == "LUMP", ]

#### MNG 
MNG.me <- data.frame(bval.min = apply(bval.MNG, 1, min, na.rm=T),
                     MNG.me = NA,
                     cg = rownames(bval.MNG), 
                     stringsAsFactors = F)
MNG.me$MNG.me <- ifelse(MNG.me$bval.min > 0.3, "LUMP", "No")
MNG.me <- MNG.me[MNG.me$MNG.me == "LUMP", ]



####   LUMP
load( file = "./GSE35069_leukcyte/bval.RData")

bval_GSE35069 <- bval_GSE35069[rownames(bval),]

LUMP <- data.frame(bval.min = apply(bval_GSE35069, 1, max, na.rm=T),
                   LUMP = NA, 
                   cg = rownames(bval_GSE35069), 
                   stringsAsFactors = F)
LUMP$LUMP <- ifelse(LUMP$bval.min <0.05, "LUMP", "No")
LUMP <- LUMP[LUMP$LUMP == "LUMP", ]

###LeMoN score 
LeMoN.cg  <- intersect(intersect(CNS.me$cg, LUMP$cg), MNG.me$cg)

load("./GSE90496_Capper/bval_MNG.RData")

bval.MNG.lemon <- bval.MNG[LeMoN.cg , ]
colMeans(bval.MNG.lemon, na.rm = T)


####heatmap 
library(pheatmap)
bval.CNS.lemon <- bval[LeMoN.cg,]
bval.heatmap <- cbind(bval.CNS.lemon, bval_GSE35069[LeMoN.cg,])

annotation_col = data.frame(
  type = c(rep("CNS",2801 ), rep("Leukocyte",60) )
)
rownames(annotation_col) = colnames(bval.heatmap)

pdf("./plots/102Heatmap_LeMoN_CNS.pdf", width = 7, height = 5)
#breaksList = seq(-3, 3, by = 0.01)
pheatmap(bval.heatmap, 
         #scale = "row",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = F, show_colnames = F,
         fontsize_row = 1, fontsize_col = 1,
         cluster_cols = T,
         cluster_rows = T,
         border_color = NA,
         color = colorRampPalette(c("darkblue", "black", "yellow"))(50), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         #breaks = breaksList,
         annotation_legend = T,
         annotation_col = annotation_col,
         annotation_colors  = list(
           type = c(CNS =pal.type[1], Leukocyte =pal.type[2])
         ),
         main = "LeMoN Score"
)
dev.off()

##MNG
bval.heatmap <- cbind(bval.MNG.lemon, bval_GSE35069[LeMoN.cg,])

annotation_col = data.frame(
  Grade = c(sampleinfo.MNG$WHO.Grade, rep("Leukocyte",60) )
)
rownames(annotation_col) = colnames(bval.heatmap)

pdf("./plots/102Heatmap_LeMoN_MNG.pdf", width = 7, height = 5)
#breaksList = seq(-3, 3, by = 0.01)
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
           Grade = c(I =pal.grade[1], II =pal.grade[2], III = pal.grade[3],
                     Leukocyte = "#808080")
         ),
         main = "LeMoN Score"
)
dev.off()



#correlation of LeMoN score with Purity
sampleinfo.MNG$Estimated.tumour.purity.according.to.TCGA...Ceccarelli.et.al..2016
LeMoN.score <- data.frame(lemon.score = colMeans(bval.MNG.lemon, na.rm = T),
                          purity = sampleinfo.MNG$Estimated.tumour.purity.according.to.TCGA...Ceccarelli.et.al..2016,
                          stringsAsFactors = F)

library(grid)
library(ggplot2)
grob3 = grobTree(textGrob(paste("Pearson Correlation:", 
                                round(cor(LeMoN.score$lemon.score, LeMoN.score$purity), 4) ), 
                          x =0.5, y = 0.1, hjust = 0, 
                          gp = gpar(col = "red", fontsize = 10, fontface = "bold")))
grob4 = grobTree(textGrob(paste("p = ", round(cor.test(LeMoN.score$lemon.score, LeMoN.score$purity)$p.value, 4) ), 
                          x =0.5, y = 0.05, hjust = 0, 
                          gp = gpar(col = "red", fontsize = 10, fontface = "bold")))

pdf("./plots/102Correlation_LeMoNpurity.pdf", width = 4, height = 4.5)
ggplot(LeMoN.score, aes(x=lemon.score, y=purity)) + 
  geom_point() + 
  ggtitle("LeMoN score vs Purity") + 
  geom_smooth(method=lm, se=FALSE) + 
  scale_x_continuous(name = "LeMoN score", 
                     limits = c(0.4, 1), 
                     breaks = seq(0, 1, 0.2)) + 
  scale_y_continuous(name = "Purity", 
                     limits = c(0.2, 0.8), 
                     breaks = seq(0, 1, 0.2)) + 
  annotation_custom(grob3) + annotation_custom(grob4) + 
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_blank(), 
        axis.line = element_line(color="black"), 
        axis.line.x = element_line(color="black"))
dev.off()
