rm(list = ls())

#library(limma)
library(RColorBrewer)
pal.grade <- brewer.pal(8, "Dark2")
pal.batch <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Set1")

load("./02eSet_subtype.RData")
load(file = "./sampleinfo_RNA.RData")
load(file = "./03DEG_subtype.RData")

subtype3.gene <- DEG[DEG$subtype == "subtype3",]$DEG
eSet <- as.data.frame(eSet)
eSet.sub3 <- eSet[subtype3.gene,]

RF.mat <- as.data.frame(t(eSet.sub3))
RF.mat$subtype <- sampleinfo$subtype

#RF.mat$subtype <- ifelse(RF.mat$subtype == "subtype3", 1, 0)
RF.mat$subtype <- as.numeric( substr(RF.mat$subtype, start = 8L, stop = 8L) )

######Random Forest######

library("randomForest")
print("---PERFORMING RANDOM FOREST ---")
set.seed(12321)
RF.output <- randomForest(subtype ~ ., data = RF.mat, importance = T, 
                          na.action = na.omit, 
                          proximity = T)
print("---TEST DONE---")
print(RF.output)
print(importance(RF.output, type=1))
summary(RF.output)

pdf("./plots/07RFvarImp.pdf", width = 7, height = 7)
varImpPlot(RF.output, main = "Importance of variables", 
           cex = 1)
dev.off()

set.seed(1234)
RF.cv <- replicate(5, rfcv(RF.mat[-ncol(RF.mat)], RF.mat$subtype, 
                           cv.fold = 10, step = 1.5), 
                   simplify = F)
RF.cv 
RF.cv <- data.frame(sapply(RF.cv, '[[', 'error.cv'))
RF.cv$gene <- rownames(RF.cv)
RF.cv <- reshape2::melt(RF.cv, id = 'gene')
RF.cv$gene <- as.numeric(as.character(RF.cv$gene))

library(ggplot2)
library(splines)
pdf("./plots/07RF_cv.pdf", width = 6, height = 5)
ggplot(RF.cv, aes(gene, value))+
  geom_smooth(se = F, method = "glm", formula = y~ns(x,6))+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color='black', 
                                        fill = 'transparent'),
        axis.text.y = element_text(face="bold", 
                                   #color="#993333", 
                                   size=14, angle=0),
        axis.text.x = element_text(face="bold", 
                                   #color="#993333", 
                                   size=12, angle=70),
        axis.text=element_text(size=14))+
  labs(title = '', 
       x = 'Numbers of Genes', 
       y = 'Cross-validation error', cex = 14)+
  scale_x_continuous(breaks=seq(0, 55, 2))
dev.off()

#ALPL
#TIMP3
#INMT
#SLC16A9

save(RF.output, RF.cv, file = "./07RF.RData" )


######ROC ######
library(pROC)
RF.mat$recurrence <- RF.mat$subtype
RF.mat$recurrence <- gsub(1, 0, RF.mat$recurrence)
RF.mat$recurrence <- gsub(2, 0, RF.mat$recurrence)
RF.mat$recurrence <- gsub(4, 0, RF.mat$recurrence)
RF.mat$recurrence <- gsub(3, 1, RF.mat$recurrence)
RF.mat$recurrence <- as.numeric(RF.mat$recurrence)

pdf("./plots/07ROC_gene_subtype.pdf", width = 4, height = 4)
modelroc <- roc(RF.mat$recurrence, RF.mat$ALPL)
plot(modelroc, print.auc=TRUE, auc.polygon=F, grid=c(0.2, 0.2),
     #grid.col=c("green", "red"),
     max.auc.polygon=F,
     auc.polygon.col="skyblue", 
     print.thres=F,
     xlim=c(1,0),
     main = "ROC for prediction of subtype by ALPL")

modelroc <- roc(RF.mat$recurrence, RF.mat$TIMP3)
plot(modelroc, print.auc=TRUE, auc.polygon=F, grid=c(0.2, 0.2),
     #grid.col=c("green", "red"),
     max.auc.polygon=F,
     auc.polygon.col="skyblue", 
     print.thres=F,
     main = "ROC for prediction of subtype by TIMP3")

modelroc <- roc(RF.mat$recurrence, RF.mat$INMT)
plot(modelroc, print.auc=TRUE, auc.polygon=F, grid=c(0.2, 0.2),
     #grid.col=c("green", "red"),
     max.auc.polygon=F,
     auc.polygon.col="skyblue", 
     print.thres=F,
     main = "ROC for prediction of subtype by INMT")

modelroc <- roc(RF.mat$recurrence, RF.mat$SLC16A9)
plot(modelroc, print.auc=TRUE, auc.polygon=F, grid=c(0.2, 0.2),
     #grid.col=c("green", "red"),
     max.auc.polygon=F,
     auc.polygon.col="skyblue", 
     print.thres=F,
     main = "ROC for prediction of subtype by SLC16A9")

dev.off()



######Violin plot ######
rownames(eSet.sub3)
feature.gene <- c("ALPL", "TIMP3", "INMT", "SLC16A9")
violin.df <- as.data.frame(t(eSet.sub3[feature.gene,]))
violin.df$subtype <- sampleinfo$subtype
colnames(violin.df)


aov(violin.df$ALPL~violin.df$subtype)
summary(aov(violin.df$ALPL~violin.df$subtype))
TukeyHSD(aov(violin.df$ALPL~violin.df$subtype))

library(ggplot2)
library(tidyverse)
pdf("./plots/07Violinplot_RF.gene.pdf", width = 5, height = 5)
violin.df %>%
  ggplot(aes(x= subtype, y =ALPL, fill = subtype))+
  geom_boxplot(outlier.size = 0)+
  geom_jitter(width = 0.3, alpha = 0.7)+
  xlab("Subtype")+
  scale_fill_manual(values=pal.subtype)+
  scale_color_manual(values=pal.subtype)+
  #facet_wrap(~miRID)+
  #theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme_classic()+
  ylab("log2(TPM)")+
  ggtitle("ALPL")

violin.df %>%
  ggplot(aes(x= subtype, y =TIMP3, fill = subtype))+
  geom_boxplot(outlier.size = 0)+
  geom_jitter(width = 0.3, alpha = 0.7)+
  xlab("Subtype")+
  scale_fill_manual(values=pal.subtype)+
  scale_color_manual(values=pal.subtype)+
  #facet_wrap(~miRID)+
  #theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme_classic()+
  ylab("log2(TPM)")+
  ggtitle("TIMP3")
violin.df %>%
  ggplot(aes(x= subtype, y =INMT, fill = subtype))+
  geom_boxplot(outlier.size = 0)+
  geom_jitter(width = 0.3, alpha = 0.7)+
  xlab("Subtype")+
  scale_fill_manual(values=pal.subtype)+
  scale_color_manual(values=pal.subtype)+
  #facet_wrap(~miRID)+
  #theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme_classic()+
  ylab("log2(TPM)")+
  ggtitle("INMT")
violin.df %>%
  ggplot(aes(x= subtype, y =SLC16A9, fill = subtype))+
  geom_boxplot(outlier.size = 0)+
  geom_jitter(width = 0.3, alpha = 0.7)+
  xlab("Subtype")+
  scale_fill_manual(values=pal.subtype)+
  scale_color_manual(values=pal.subtype)+
  #facet_wrap(~miRID)+
  #theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme_classic()+
  ylab("log2(TPM)")+
  ggtitle("SLC16A9")
dev.off()
