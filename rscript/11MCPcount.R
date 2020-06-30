rm(list = ls())

library(RColorBrewer)
pal.grade <- brewer.pal(8, "Dark2")
pal.batch <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Set1")

#load(file = "./TPM.RData")
load("./01TPM_gene.RData")
load("./02eSet_subtype.RData")
load(file = "./sampleinfo_RNA.RData")

library(MCPcounter)
ExampleEstimates=MCPcounter.estimate(MCPcounterExampleData,
                                     featuresType="affy133P2_probesets")

MCP.out <- MCPcounter.estimate(eSet,featuresType="HUGO_symbols")
MCP.out <- as.matrix(MCP.out)
rownames(MCP.out)
sampleinfo <- cbind(sampleinfo, as.data.frame(t(MCP.out)))

### heatmap
library(pheatmap)

annotation_col = data.frame(
  Grade = sampleinfo$grade,
  #Batch = sampleinfo$batch,
  Subtype = sampleinfo$subtype
)
rownames(annotation_col) = colnames(MCP.out)
breaksList = seq(-1, 1, by = 0.1)
pdf("./plots/11Heatmap_MCP.pdf", width = 5.5, height = 4)
pheatmap(MCP.out[1:9,order(sampleinfo$subtype)], 
         scale = "column",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = T, show_colnames = F,
         fontsize_row = 12, fontsize_col = 1,
         cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         #breaks = breaksList,
         annotation_legend = T,
         annotation_col = annotation_col,
         annotation_colors  = list(
           #Batch = c(GSE85133 =pal.batch[1], GSE136661 =pal.batch[2]),
           Subtype = c(subtype1 = pal.subtype[1], subtype2=pal.subtype[2], 
                       subtype3 = pal.subtype[3], subtype4 = pal.subtype[4]),
           Grade = c(I = pal.grade[1], II =pal.grade[2], III=pal.grade[3])),
         main = "Immune Cells Infiltration \nin Meningioma Subtypes"
)
dev.off()




#### MCP with mRNAsi 



####  MCP between subtypes 
library(ggplot2)

violin.df <- data.frame(Grade = sampleinfo$grade, 
                        Subtype = sampleinfo$subtype, 
                        StemIndex = sampleinfo$mRNAsi, 
                        Gender = sampleinfo$gender,
                        Age_number = sampleinfo$age,
                        stringsAsFactors = F)
violin.df <- cbind(violin.df, t(MCP.out))
colnames(violin.df) <- gsub(" ","_", colnames(violin.df))


pdf("./plots/11Violinplot_immuneCell.pdf", width = 5, height = 5)
summary(aov(Cytotoxic_lymphocytes ~ Subtype, data = violin.df))
TukeyHSD(aov(Cytotoxic_lymphocytes ~ Subtype, data = violin.df))

ggplot(data = violin.df, aes(x = Subtype, y = Cytotoxic_lymphocytes))+
  geom_violin(aes(fill = factor(Subtype)))+
  geom_boxplot(width = 0.3, outlier.size = -1)+
  geom_jitter(height = 0.002, width = 0.2, cex = 1.2)+
  scale_fill_brewer(palette="Set1")+
  labs(y = "Cytotoxic Lymphocytes (%)", 
       x = "ANOVA: p = 3.76e-12 ", 
       title = "")+ 
  #ylim(0,1)+
  #theme_bw()+ 
  theme_classic()+
  theme(legend.position = "none") 


summary(aov(NK_cells ~ Subtype, data = violin.df))
TukeyHSD(aov(NK_cells ~ Subtype, data = violin.df))

ggplot(data = violin.df, aes(x = Subtype, y = NK_cells))+
  geom_violin(aes(fill = factor(Subtype)))+
  geom_boxplot(width = 0.3, outlier.size = -1)+
  geom_jitter(height = 0.002, width = 0.2, cex = 1.2)+
  scale_fill_brewer(palette="Set1")+
  labs(y = "NK cells (%) ", 
       x = "ANOVA: p = 0.000984", 
       title = "")+ 
  #ylim(0,1)+
  #theme_bw()+ 
  theme_classic()+
  theme(legend.position = "none") 


summary(aov(Monocytic_lineage ~ Subtype, data = violin.df))
TukeyHSD(aov(Monocytic_lineage ~ Subtype, data = violin.df))

ggplot(data = violin.df, aes(x = Subtype, y = Monocytic_lineage))+
  geom_violin(aes(fill = factor(Subtype)))+
  geom_boxplot(width = 0.3, outlier.size = -1)+
  geom_jitter(height = 0.002, width = 0.2, cex = 1.2)+
  scale_fill_brewer(palette="Set1")+
  labs(y = "Monocytic Lineage (%) ", 
       x = "ANOVA: p = 1.6e-05 ", 
       title = "")+ 
  #ylim(0,1)+
  #theme_bw()+ 
  theme_classic()+
  theme(legend.position = "none") 

dev.off()



######PDL1#######
eSet <- log2(TPM_gene+1)
grep("CD274", rownames(eSet))
PDL1.df <- as.data.frame(eSet["CD274",])
PDL1.df <- cbind(PDL1.df, sampleinfo$subtype)
colnames(PDL1.df)[1] <- "PDL1"
colnames(PDL1.df)[2] <- "subtype"

library(ggplot2)

summary(aov(PDL1 ~ subtype, data = PDL1.df))
TukeyHSD(aov(PDL1 ~ subtype, data = PDL1.df))

pdf("./plots/11Violinplot_PDL1.pdf", width = 4, height = 4)
ggplot(data = PDL1.df, aes(x = subtype, y = PDL1))+
  geom_violin(aes(fill = factor(subtype)))+
  geom_boxplot(width = 0.3, outlier.size = -1)+
  geom_jitter(height = 0.002, width = 0.2, cex = 1.2)+
  scale_fill_brewer(palette="Set1")+
  labs(y = "PDL1 (log2 TPM)", 
       x = "ANOVA: p = 0.000891", 
       title = "")+ 
  #ylim(0,1)+
  #theme_bw()+ 
  theme_classic()+
  theme(legend.position = "none") 
dev.off()
