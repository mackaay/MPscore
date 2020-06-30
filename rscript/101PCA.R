rm(list = ls())

load("./GSE90496_Capper/bval_MNG.RData")
library(limma)
library(RColorBrewer)
sampleinfo.MNG$WHO.Grade
sampleinfo.MNG$Pathological.Diagnosis..WHO.2016.


library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
keep <- !(rownames(bval.MNG) %in% ann450k$Name[ann450k$chr %in% 
                                                      c("chrX","chrY")])


pal.grade <- brewer.pal(8, "Dark2")


pdf("./plots/01PCA_GSE90946.pdf", width = 5, height = 5)
plotMDS(bval.MNG, top=1000, gene.selection="common", 
        col=pal.grade[factor(sampleinfo.MNG$WHO.Grade)], 
        pch = 19, cex = 1.2,
        main = "WHO Grade")
legend("bottomleft", legend=levels(factor(sampleinfo.MNG$WHO.Grade)), 
       text.col=pal.grade,
       bg="white", cex=1)

plotMDS(bval.MNG, top=1000, gene.selection="common", 
        col=pal.grade[factor(sampleinfo.MNG$WHO.Grade)], 
        pch = 19, cex = 1.2, dim= c(1,3),
        main = "WHO Grade")
legend("topleft", legend=levels(factor(sampleinfo.MNG$WHO.Grade)), 
       text.col=pal.grade,
       bg="white", cex=1)

plotMDS(bval.MNG, top=1000, gene.selection="common", 
        col=pal.grade[factor(sampleinfo.MNG$WHO.Grade)], 
        pch = 19, cex = 1.2, dim= c(2,3),
        main = "WHO Grade")
legend("topleft", legend=levels(factor(sampleinfo.MNG$WHO.Grade)), 
       text.col=pal.grade,
       bg="white", cex=1)
dev.off()

