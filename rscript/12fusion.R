rm(list = ls())
options(stringsasfactors = F)

###GSE85133
dir <- "/datasets/work/hb-diab-cfdna/work/scratch/chenkai/AnaM/GSE85133_RNAseq_raw/alignment_fusion/"
fusion.files <- list.files(path = dir,pattern = '*predictions.tsv$', recursive = T)


fusion <- read.table(paste(dir, fusion.files[1], sep = ""), 
                     sep = '\t', header = F, stringsAsFactors = F)
fusion <- fusion[,c(1,2,5:8,11,12,17)]
colnames(fusion) <- c("FusionName", "ReadCount", "LeftGene", "LeftBreakpoint", 
                      "RightGene", "RightBreakpoint", "LargeAnchorSupport", "FFPM", "anno")
fusion$sample <- rep(substr(fusion.files[1], start = 1L, stop = 10L),nrow(fusion))

for (i in 2:length(fusion.files)) {
  tmp <- read.table(paste(dir, fusion.files[i], sep = ""), 
                       sep = '\t', header = F, stringsAsFactors = F)
  tmp <- tmp[,c(1,2,5:8,11,12,17)]
  colnames(tmp) <- c("FusionName", "ReadCount", "LeftGene", "LeftBreakpoint", 
                        "RightGene", "RightBreakpoint", "LargeAnchorSupport", "FFPM", "anno")
  tmp$sample <- rep(substr(fusion.files[i], start = 1L, stop = 10L),nrow(tmp))
  
  fusion <- rbind(fusion, tmp)
}
table(fusion$sample)
table(fusion$FusionName)
table(c(fusion$LeftGene, fusion$RightGene))



###GSE136661
dir <- "/datasets/work/hb-diab-cfdna/work/scratch/chenkai/AnaM/GSE136661_RNAseq_raw/alignment_fusion/"
fusion.files <- list.files(path = dir,pattern = '*predictions.tsv$', recursive = T)

fusion2 <- read.table(paste(dir, fusion.files[1], sep = ""), 
                     sep = '\t', header = F, stringsAsFactors = F)
fusion2 <- fusion2[,c(1,2,5:8,11,12,17)]
colnames(fusion2) <- c("FusionName", "ReadCount", "LeftGene", "LeftBreakpoint", 
                      "RightGene", "RightBreakpoint", "LargeAnchorSupport", "FFPM", "anno")
fusion2$sample <- rep(substr(fusion.files[1], start = 1L, stop = 11L),nrow(fusion2))

for (i in c(2:14, 16:31, 33:60, 62:length(fusion.files))) {
  tmp <- read.table(paste(dir, fusion.files[i], sep = ""), 
                    sep = '\t', header = F, stringsAsFactors = F)
  tmp <- tmp[,c(1,2,5:8,11,12,17)]
  colnames(tmp) <- c("FusionName", "ReadCount", "LeftGene", "LeftBreakpoint", 
                     "RightGene", "RightBreakpoint", "LargeAnchorSupport", "FFPM", "anno")
  tmp$sample <- rep(substr(fusion.files[i], start = 1L, stop = 11L),nrow(tmp))
  
  fusion2 <- rbind(fusion2, tmp)
}
table(fusion2$sample)
table(fusion2$FusionName)
sort(table(c(fusion2$LeftGene, fusion2$RightGene)))


#merge two dataset together
fusion_merge <- rbind(fusion, fusion2)

write.csv(fusion_merge, file = "./fusion.csv",sep = ",", row.names = F)
#fusion_merge <- read.csv("./fusion.csv", stringsAsFactors = F)

table(fusion_merge$sample)
fusion.df <- as.data.frame(table(fusion_merge$sample), stringsAsFactors = F)
nofusion <- data.frame(Var1 = c("SRR10042635","SRR10042652", "SRR10042681" ), 
                       Freq = rep(0, 3))
fusion.df <- rbind(fusion.df, nofusion)
fusion.df <- fusion.df[order(fusion.df$Var1),]
fusion.df <- fusion.df[c(161:179, 1:160),]

load(file = "./sampleinfo_RNA.RData")
sampleinfo$fusion.freq <- fusion.df$Freq
save(sampleinfo, file = "./sampleinfo_RNA.RData")

library(RColorBrewer)
pal.grade <- brewer.pal(8, "Dark2")
pal.batch <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Set1")

pdf("./plots/12Boxplot_fusionFreq.pdf", width = 4, height = 5)
boxplot(sampleinfo$fusion.freq ~ sampleinfo$subtype, col = pal.subtype, 
        main = "Gene Fusion Frequency", 
        xlab = "Subtypes", ylab = "Number of Fusions")
boxplot(sampleinfo$fusion.freq ~ sampleinfo$grade, col = pal.grade, 
        main = "Gene Fusion Frequency", 
        xlab = "Grade", ylab = "Number of Fusions")
dev.off()

#subtypes ANOVA
summary(aov(sampleinfo$fusion.freq ~ sampleinfo$subtype))
TukeyHSD(aov(sampleinfo$fusion.freq ~ sampleinfo$subtype))
#grade ANOVA
summary(aov(sampleinfo$fusion.freq ~ sampleinfo$grade))
TukeyHSD(aov(sampleinfo$fusion.freq ~ sampleinfo$grade))



###fusion feature in each subtype
sort(table(fusion_merge$FusionName))

sampleinfo[sampleinfo$sra %in% c("SRR10042635","SRR10042652", "SRR10042681" ),]


#
subtype.sra <- sampleinfo[sampleinfo$subtype == "subtype3",]$sra 
fusion.subtype <- fusion_merge[fusion_merge$sample %in% subtype.sra,]$FusionName
sort(table(fusion.subtype))


#NF
fusion.NF <- as.data.frame(sort(table(fusion_merge$FusionName)), stringsAsFactors = F)
grep("NF2", fusion.NF$Var1, value = T)


