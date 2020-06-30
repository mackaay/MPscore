rm(list = ls())

#library(limma)
library(RColorBrewer)
pal.grade <- brewer.pal(8, "Dark2")
pal.batch <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Set1")


######GSE74385 recurrence vs nonrecurrence######
load(file = "./00GSE74385_validation.RData")

eSet["ALPL",]

violin.df <- data.frame(Type = sampleinfo$`subtype:ch1`, 
                        ALPL = eSet["ALPL",], 
                        stringsAsFactors = F)

library(ggplot2)
pdf("./plots/09Violinplot_subtype3gene.pdf", width = 5, height = 5)
summary(aov(ALPL ~ Type, data = violin.df))
TukeyHSD(aov(ALPL ~ Type, data = violin.df))

ggplot(data = violin.df, aes(x = Type, y = ALPL))+
  geom_violin(aes(fill = factor(Type)))+
  geom_boxplot(width = 0.3, outlier.size = -1)+
  geom_jitter(height = 0.002, width = 0.2, cex = 1.2)+
  scale_fill_brewer(palette="Dark2")+
  labs(y = "Stem Index", 
       x = "ANOVA: p = 0.00195 ", 
       title = "")+ 
  #ylim(0,1)+
  #theme_bw()+ 
  theme_classic()+
  theme(legend.position = "none") 


######GSE16581  survival #######
rm(list = ls())
#library(limma)
library(RColorBrewer)
pal.grade <- brewer.pal(8, "Dark2")
pal.batch <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Set1")

load(file = "./00GSE16581_array.RData")

violin.df <- as.data.frame(t(eSet[eSet$symbol == "ALPL",]))
violin.df <- violin.df[1:68,]
violin.df$`1557924_s_at` <- as.numeric(as.character(violin.df$`1557924_s_at`))
violin.df$`215783_s_at` <- as.numeric(as.character(violin.df$`215783_s_at`))
violin.df$Recurrence <- sampleinfo$`recurrence_frequency:ch1`
violin.df$Recurrence <- ifelse(violin.df$Recurrence == 0 , "no", "yes")
colnames(violin.df)[1:2] <- c("ALPL_1", "ALPL_2")
violin.df$grade <- sampleinfo$`who grade:ch1`
violin.df$Rfreq <- sampleinfo$`recurrence_frequency:ch1`
violin.df$Rfreq <- as.numeric(violin.df$Rfreq)
violin.df$ALPL <- rowMeans(cbind(log2(violin.df$ALPL_1), log2(violin.df$ALPL_2)))


###   Correlation
cor(violin.df$ALPL, violin.df$Rfreq)
cor.test(violin.df$ALPL, violin.df$Rfreq)
library(grid)
SPP1 <- violin.df$ALPL
CALB2 <- violin.df$Rfreq
grob3 = grobTree(textGrob(paste("Pearson Correlation:", 
                                round(cor(CALB2, SPP1), 4) ), 
                          x =0.5, y = 0.1, hjust = 0, 
                          gp = gpar(col = "red", fontsize = 10, fontface = "bold")))
grob4 = grobTree(textGrob(paste("p = ", round(cor.test(CALB2, SPP1)$p.value, 4) ), 
                          x =0.5, y = 0.05, hjust = 0, 
                          gp = gpar(col = "red", fontsize = 10, fontface = "bold")))
pdf("./plots/09Cor_ALPLvRfeq", width = 5, height = 5)
ggplot(violin.df, aes(x=ALPL, y=Rfreq)) + 
  geom_point() + 
  ggtitle("ALPL vs Recurrence Frequncy") + 
  geom_smooth(method=lm, se=FALSE) + 
  scale_x_continuous(name = "ALPL", 
                     limits = c(3, 10), 
                     breaks = seq(3, 10, 1)) + 
  scale_y_continuous(name = "Recurrence Frequncy", 
                     limits = c(0, 10), 
                     breaks = seq(0, 10, 1)) + 
  annotation_custom(grob3) + 
  annotation_custom(grob4) + 
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_blank(), 
        axis.line = element_line(color="black"), 
        axis.line.x = element_line(color="black"))
dev.off()

### Violin plot
library(ggplot2)
#summary(aov(log2(ALPL_1) ~ Recurrence, data = violin.df))
#TukeyHSD(aov(ALPL_1 ~ Recurrence, data = violin.df))

pdf("./plots/09Violinplot_ALPL.pdf", width = 5, height = 5)

t.test(violin.df$ALPL ~ violin.df$Recurrence)
ggplot(data = violin.df, aes(x = Recurrence, y = ALPL))+
  geom_violin(aes(fill = factor(Recurrence)))+
  geom_boxplot(width = 0.3, outlier.size = -1)+
  geom_jitter(height = 0.002, width = 0.2, cex = 1.2)+
  scale_fill_brewer(palette="Dark2")+
  facet_grid(~grade)+
  labs(y = "ALPL", 
       x = "t test: p = 0.0001193 ", 
       title = "ALPL expression by WHO Grade")+ 
  #ylim(0,1)+
  #theme_bw()+ 
  theme_classic()+
  theme(legend.position = "none") 
dev.off()



####   Survival
violin.df$status <- ifelse(sampleinfo$`vital status:ch1` == "Alive", 0, 1)
violin.df$time <- as.numeric(sampleinfo$`tts:ch1`)
violin.df$group <- ifelse(violin.df$ALPL > mean(violin.df$ALPL), "lowRisk", "highRisk")


library(survival)
my.surv <- Surv(violin.df$time,violin.df$status==1)
kmfit1 <- survfit(my.surv~group,data = violin.df) 

pdf("./plots/107KMplot_PCR.pdf", width = 4, height = 3)
plot(kmfit1,col = c("#ffdb0d","#0daeff"), 
     main = "Overall Survival", 
     xlab = "Days", 
     ylab = "Cumulative Survival Rate")
legend("bottomleft", legend = c("LRRC39 high", "LRRC39 low"), 
       col = c("#ffdb0d","#0daeff"), lty = c(1,1), cex = 0.75)
dev.off()
