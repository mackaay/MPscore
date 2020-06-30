rm(list = ls())

#library(limma)
library(RColorBrewer)
pal.grade <- brewer.pal(8, "Dark2")
pal.batch <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Set1")

load("./02eSet_subtype.RData")
load(file = "./sampleinfo_RNA.RData")
