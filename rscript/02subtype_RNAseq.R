rm(list = ls())

library(limma)
library(RColorBrewer)

load("./01TPM_gene.RData")

eSet <- removeBatchEffect(log2(TPM_gene+1), sampleinfo$batch)


######subtype clustering######
library(ConsensusClusterPlus)

eSet_cluster = sweep(eSet,1, apply(eSet,1,median,na.rm=T)) #median center dataset
#mVals_5000 <- na.omit(mVals_5000)

title= c("./Subtype_cluster/")
eSet_cluster <- as.matrix(eSet_cluster)
results = ConsensusClusterPlus(eSet_cluster, maxK=10,reps=50,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",
                               seed=71279,plot="pdf")
results[[4]][["consensusMatrix"]][1:5,1:5]
results[[4]][["consensusTree"]]
results[[4]][["consensusClass"]][1:5]
consensus.Mat <- results[[4]][["consensusMatrix"]]
save(consensus.Mat ,  file = "./Subtype_cluster/02ConsensusClusterPlus_Matrix.RData")

icl = calcICL(results,title=title,plot="pdf")
icl[["clusterConsensus"]]
icl[["itemConsensus"]][1:5,]

#clustering k = 4 is ideal subgroup for downstream analysis
AnaM_subtype <- results[[4]][["consensusClass"]]

sampleinfo$subtype <- AnaM_subtype


pal.grade <- brewer.pal(8, "Dark2")
pal.batch <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Set1")


pdf("./plots/02PCA_RNAsubtype.pdf", width = 5, height = 5)
plotMDS(eSet, top=1000, gene.selection="common", 
        col=pal.subtype[factor(sampleinfo$subtype)], 
        pch = 19, cex = 1.2,
        main = "Subtype")
legend("bottomleft", legend=levels(factor(sampleinfo$subtype)), text.col=pal.subtype,
       bg="white", cex=1)
dev.off()


save(eSet, sampleinfo, file = "./02eSet_subtype.RData")



######SigClust Example######
## Simulate a dataset from a collection of mixtures of two
## multivariate Gaussian distribution with different means.
library(sigclust)
mu <- 5
n <- 30
p <- 500
dat <- matrix(rnorm(p*2*n),2*n,p)
dat[1:n,1] <- dat[1:n,1]+mu
dat[(n+1):(2*n),1] <- dat[(n+1):(2*n),1]-mu

nsim <- 1000
nrep <- 1
icovest <- 3
pvalue <- sigclust(dat,nsim=nsim,nrep=nrep,labflag=0,icovest=icovest)
#sigclust plot
plot(pvalue)

######CancerSubtype SigClust######
library(CancerSubtypes)
#data(GeneExp)
#data(miRNAExp)
#data(time)
#data(status)
#GBM=list(GeneExp=GeneExp,miRNAExp=miRNAExp)
#result=ExecuteSNF(GBM, clusterNum=3, K=20, alpha=0.5, t=20)
#group=result$group
#sigclust1=sigclustTest(miRNAExp,group, nsim=500, nrep=1, icovest=3)
#sigclust2=sigclustTest(miRNAExp,group, nsim=1000, nrep=1, icovest=1)

group <- sampleinfo$subtype
sigclust1 <- sigclustTest(eSet,group, nsim=500, nrep=1, icovest=3)
sigclust2 <- sigclustTest(eSet_cluster,group, nsim=500, nrep=1, icovest=3)

#save(group, eSet, eSet_cluster, file = "./Subtype_cluster/02CancerSubtype_sigclust.RData")



######Silhouette width######
library(CancerSubtypes)
library(RColorBrewer)
load("./Subtype_cluster/02CancerSubtype_sigclust.RData")
load("./Subtype_cluster/02ConsensusClusterPlus_Matrix.RData")
#data(GeneExp)
#data(miRNAExp)
#GBM=list(GeneExp=GeneExp,miRNAExp=miRNAExp)
#result=ExecuteSNF(GBM, clusterNum=3, K=20, alpha=0.5, t=20,plot = FALSE)
group.s <- group[order(group)]
consensus.Mat.s <- consensus.Mat[order(group),order(group)]
sil=silhouette_SimilarityMatrix(group.s, consensus.Mat.s)
pdf("./plots/02SilhouettePlot.pdf", width = 5, height = 4)
plot(sil , col = pal.subtype[1:4])
dev.off()



####Sankey plot######
load("./02eSet_subtype.RData")
sampleinfo$subtype <- paste("subtype", sampleinfo$subtype, sep = "")

library(ggplot2)
library(ggalluvial)
sankey.df <- sampleinfo[,c("grade", "subtype")]
sankey.df$freq <- rep(1, nrow(sankey.df))

pdf("./plots/02SankeyPlot.pdf", width = 5, height = 5)
ggplot(sankey.df, aes(y = freq, axis1 = grade, axis2 = subtype)) +
  geom_alluvium(aes(fill = subtype), width = 0.5) +
  geom_stratum() +
  geom_label(stat = "stratum", infer.label = T) +
  scale_x_discrete(limits = c("grade", "subtype")) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme(legend.position = "none")+
  ggtitle("Sankey plot")+
  theme_void()
dev.off()


#####Table Statistical#####
sampleinfo$age <- gsub("U", NA, sampleinfo$age)
sampleinfo$age <-  as.integer(sampleinfo$age)
summary(sampleinfo$age)

summary(sampleinfo[sampleinfo$subtype == "subtype1",])
summary(sampleinfo[sampleinfo$subtype == "subtype2",])
summary(sampleinfo[sampleinfo$subtype == "subtype3",])
summary(sampleinfo[sampleinfo$subtype == "subtype4",])

sampleinfo$gender
table(sampleinfo[sampleinfo$subtype == "subtype1",]$gender)
table(sampleinfo[sampleinfo$subtype == "subtype2",]$gender)
table(sampleinfo[sampleinfo$subtype == "subtype3",]$gender)
table(sampleinfo[sampleinfo$subtype == "subtype4",]$gender)

table(sampleinfo[sampleinfo$subtype == "subtype1",]$grade)
table(sampleinfo[sampleinfo$subtype == "subtype2",]$grade)
table(sampleinfo[sampleinfo$subtype == "subtype3",]$grade)
table(sampleinfo[sampleinfo$subtype == "subtype4",]$grade)


