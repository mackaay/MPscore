rm(list = ls())

library(limma)
library(RColorBrewer)
pal.grade <- brewer.pal(8, "Dark2")
pal.batch <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Set1")

load( file = "./03DEG_subtype.RData")
load(file = "./03DEGup_subtype.RData")

DEGup.list <- DEGup$DEG
DEG.list <- DEG$DEG



load("./00GSE84263_array.RData")

eSet.valid <- eSet[which(rownames(eSet) %in% DEG$DEG),]
eSet.valid <- eSet.valid[match(DEG$DEG, rownames(eSet.valid)),]
rownames(eSet.valid)

library(pheatmap)

annotation_col = data.frame(
  Grade = sampleinfo$grade,
  Batch = sampleinfo$batch,
  Subtype = sampleinfo$subtype
)
rownames(annotation_col) = colnames(eSet)

pdf("./plots/03Heatmap_arraySubtype.pdf", width = 7, height = 7)
breaksList = seq(-4, 4, by = 0.01)
pheatmap(eSet.valid, 
         scale = "row",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = F, show_colnames = F,
         fontsize_row = 1, fontsize_col = 1,
         cluster_cols = T,
         cluster_rows = F,
         color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         #annotation_legend = T,
         #annotation_col = annotation_col,
         #annotation_colors  = list(
           #Batch = c(GSE85133 =pal.batch[1], GSE136661 =pal.batch[2]),
           #Subtype = c(subtype1 = pal.subtype[1], subtype2=pal.subtype[2], 
            #           subtype3 = pal.subtype[3], subtype4 = pal.subtype[4]),
           #Grade = c(I = pal.grade[1], II =pal.grade[2], III=pal.grade[3])),
         main = "GSE84263"
)
dev.off()



######GSE74385 array ######
load("./00GSE74385_array.RData")
eSet.DEGup <- eSet_GSE74385[match(DEGup.list, rownames(eSet_GSE74385)), ]
eSet.DEGup <- eSet.DEGup[grep("NA", rownames(eSet.DEGup), invert = T),]

eSet.DEG <- eSet_GSE74385[match(DEG.list, rownames(eSet_GSE74385)), ]
eSet.DEG <- eSet.DEG[grep("NA", rownames(eSet.DEG), invert = T),]

library(pheatmap)
col.type <- brewer.pal(9, "Paired")

annotation_col = data.frame(
  Subtype = sampleinfo_GSE74385$`subtype:ch1`
)
rownames(annotation_col) = colnames(eSet.DEGup)
names(col.type) <- levels(as.factor(sampleinfo_GSE74385$`subtype:ch1`))
breaksList = seq(-2, 2, by = 0.1)

pdf("./plots/03Heatmap_array_GSE74385.pdf",  width = 7, height = 7)
pheatmap(eSet.DEGup, 
         scale = "row",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = T, show_colnames = F,
         fontsize_row = 5, fontsize_col = 1,
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
         main = "GSE74385"
)

rownames(annotation_col) = colnames(eSet.DEG)
names(col.type) <- levels(as.factor(sampleinfo_GSE74385$`subtype:ch1`))
breaksList = seq(-2, 2, by = 0.1)
pheatmap(eSet.DEG, 
         scale = "row",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = T, show_colnames = F,
         fontsize_row = 5, fontsize_col = 1,
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
         main = "GSE74385"
)

dev.off()


library(GSVA)
gs = read.csv("./ssGSEA_geneset.csv", stringsAsFactors = FALSE, check.names = FALSE)
#a = read.table("RNA.csv", stringsAsFactors = FALSE, header = TRUE, row.names = 1, sep = ",")
a = as.matrix(eSet_GSE74385)
gs = as.list(gs)
gs = lapply(gs, function(x) x[!is.na(x)])
gs[[2]] <- gs[[2]][1:16]
gs[[3]] <- gs[[3]][1:39]

ssgsea_score = gsva(a, gs, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
ssgsea_score <- as.data.frame(t(ssgsea_score))
ssgsea_score$MPscore <- c(ssgsea_score$MP_sub3_up + (0- ssgsea_score$MP_sub3_down))
data.frame(sampleinfo_GSE74385$geo_accession, sampleinfo_GSE74385$`subtype:ch1`, 
           rownames(ssgsea_score), ssgsea_score$MPscore)
sampleinfo_GSE74385$MPscore <- ssgsea_score$MPscore
colnames(sampleinfo_GSE74385)[3] <- "type"

pdf("./plots/03Boxplot_array_MPscore.pdf", width = 4, height = 4.5)
boxplot(MPscore ~ type , data = sampleinfo_GSE74385, 
        col = col.type, notch = F, 
        ylim = c(0,1),
        lty = 1,
        las = 2, 
        main = "GSE74385")
dev.off()

t.test()

summary(aov(MPscore ~ type , data = sampleinfo_GSE74385))
TukeyHSD(aov(MPscore ~ type , data = sampleinfo_GSE74385))$type

write.csv(TukeyHSD(aov(MPscore ~ type , data = sampleinfo_GSE74385))$type, 
          file = "./03Boxplot_array_MPscore.csv", row.names = T)
