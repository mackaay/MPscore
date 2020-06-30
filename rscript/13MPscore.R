rm(list = ls())

library(RColorBrewer)
pal.grade <- brewer.pal(8, "Dark2")
pal.batch <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Set1")

#load(file = "./TPM.RData")
load("./01TPM_gene.RData")
load("./02eSet_subtype.RData")
load(file = "./sampleinfo_RNA.RData")
load( file = "./03DEG_subtype.RData")
load(file = "./03DEGup_subtype.RData")

#####MP subtype 3 score######
table(DEG$subtype)

library(GSVA)
#gs <- as.list(sample(10:100, size=100, replace=TRUE))
#gs <- lapply(gs, function(n, p) sample(1:20000, size=30, replace=FALSE), 20000) ## sample gene sets
gs = read.csv("./ssGSEA_geneset.csv", stringsAsFactors = FALSE, check.names = FALSE)
#a = read.table("RNA.csv", stringsAsFactors = FALSE, header = TRUE, row.names = 1, sep = ",")
a = as.matrix(eSet)
gs = as.list(gs)
gs = lapply(gs, function(x) x[!is.na(x)])
gs[[2]] <- gs[[2]][1:16]
gs[[3]] <- gs[[3]][1:39]

ssgsea_score = gsva(a, gs, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
ssgsea_score <- as.data.frame(t(ssgsea_score))
ssgsea_score$MPscore <- c(ssgsea_score$MP_sub3_up + (0- ssgsea_score$MP_sub3_down))
data.frame(sampleinfo$sra, sampleinfo$subtype, rownames(ssgsea_score), ssgsea_score$MPscore)
sampleinfo$MPscore <- ssgsea_score$MPscore

pdf("13Boxplot_MPscore.pdf", width = 4.5, height = 4.5)
boxplot(MPscore ~ subtype , data = sampleinfo, 
        col = pal.subtype, notch = T, 
        ylim = c(-1,1),
        main = "MP score")
dev.off()

save(sampleinfo, file = "./sampleinfo_RNA.RData")
write.csv(ssgsea_score, "ssGSEA.csv")