## DegPlotWide == modified function
library(dplyr)
library(ggpubr)

## Fig. 3: DLA vs LLA
library(DEGreport) ## to build clusters
library(ggplot2)
library(DESeq2)
load("data/DE.RData")
source("0_helper_functions.R")


## do DL and LL overlap?
intersect(row.names(DL4500.degs), row.names(LL4500.degs))
LLDLmerged <- merge(DL4500.degs, LL4500.degs, by = 0)
plot(LLDLmerged$log2FoldChange.x, LLDLmerged$log2FoldChange.y)
cor(LLDLmerged$log2FoldChange.x, LLDLmerged$log2FoldChange.y)
pLLDLmerged <- ggplot(LLDLmerged, aes(log2FoldChange.x, log2FoldChange.y)) + 
  geom_point(alpha=1) + theme_bw() + 
#  xlab("log2 fold change, 45 mM DLA") + ylab("log2 fold change, 45 mM LLA") +
  xlab(expression(Log[2]~"fold change, 45 mM DLA")) + 
  ylab(expression(Log[2]~"fold change, 45 mM LLA")) +
  geom_vline(xintercept=0, linetype="dotted") + 
  geom_hline(yintercept=0, linetype="dotted") + 
  theme(axis.title.y = element_text(hjust=.35))
pLLDLmerged

## all genes
LLDLmergedAll <- merge(res_DL45vsC.df, res_LL45vsC.df, by = 0, all = TRUE)
plot(LLDLmergedAll$log2FoldChange.x, LLDLmergedAll$log2FoldChange.y)
cor(LLDLmergedAll$log2FoldChange.x, LLDLmergedAll$log2FoldChange.y, use="complete.obs")
pLLDLmergedAll <- 
  ggplot(LLDLmergedAll, aes(log2FoldChange.x, log2FoldChange.y)) + 
  geom_point(alpha=.1) + theme_bw() + 
#  xlab("log2 fold change, 45 mM DLA") + ylab("log2 fold change, 45 mM LLA") +
  xlab(expression(Log[2]~"fold change, 45 mM DLA")) + 
  ylab(expression(Log[2]~"fold change, 45 mM LLA")) +
  geom_vline(xintercept=0, linetype="dotted") + 
  geom_hline(yintercept=0, linetype="dotted") + 
  theme(axis.title.y = element_text(hjust=.35))
pLLDLmergedAll

p3ab <- 
  ggarrange(pLLDLmerged, pLLDLmergedAll, labels = c("a", "b"))
p3ab
#ggsave("correlation_DL_LL.png", width = 20, height = 8, units = "cm")


## Fig. 3c: clusters!
degs <- 
  Reduce(union, list(row.names(DL0005.degs), row.names(DL0050.degs), row.names(DL0500.degs), 
                     row.names(DL4500.degs), row.names(LL4500.degs)))#, row.names(DLLL.degs)))
design <- as.data.frame(colData(dds))
ma = assay(rlog(dds))[row.names(res) %in% degs, ] 

#png("S2_degPatterns.png", width = 1600, height = 900, res = 200)
#res2 <- myDegPatterns(ma, design, time = "condition", minc=1, scale=F)
res2 <- degPatterns(ma, design, time = "condition", minc=1)
p3c <- res2$plot + 
  scale_x_discrete(labels = c("Control", "0.05 mM DLA", "0.5 mM DLA", "5 mM DLA", "45 mM DLA", "45 mM LLA")) + 
  theme_bw() + theme(legend.position = 'none', axis.text.x = element_text(angle = 45L, hjust = 1L))

#cluster.info <- unique(p3C$data$title)
#names(cluster.info) <- 1:6
#p3c + facet_wrap( ~ cluster, labeller = labeller(cluster=cluster.info), nrow=1) 

## Fig. S3
clust1 <- res2$df[res2$df$cluster == 1, ]
clust1$genes
writeLines(clust1$genes, "data/clust1.txt")

clust1$genes <- gsub(".A", "-A", clust1$genes)
clust1$genes <- gsub(".B", "-B", clust1$genes)

png("figs/FigS3.png", width = 20, height = 15, units="cm", res=300)
myDegPlotWide(genes=clust1$genes, counts = dds, group = "condition") + 
  expand_limits(y=0) +
  facet_wrap(~gene, nrow = 3, scales = "fixed" ) + # , labeller = labeller(gene = qual.sensor.genes)) + 
  scale_x_discrete(labels = c("Control", "0.05 mM DLA", "0.5 mM DLA", "5 mM DLA", "45 mM DLA", "45 mM LLA")) + 
  theme(strip.text = element_text(face = "italic"), legend.position = 'none') + 
  scale_color_manual(values = c("black", blues9[5:8], "yellow4"))
dev.off()

## let's now deal with the DL vs LL, trying to make sense out of these genes
degsDLvsLL <- row.names(DLLL.degs)
## those that are NOT different b/w LL and control && different from control == qual sensors
qual.sensor <- degsDLvsLL[!degsDLvsLL %in% row.names(LL4500.degs) & 
                            degsDLvsLL %in% row.names(DL4500.degs) & 
                            degsDLvsLL %in% row.names(DL0500.degs)]

qual.sensor.genes <- c("YHL028W, WSC4", "YLR121C, YPS3", "YGR189C, CRH1", "YLR054C, OSW2", "YHR209W, CRG1")
names(qual.sensor.genes) <- qual.sensor
## DLvsLL
p3d <-
  myDegPlotWide(genes=qual.sensor, counts = dds, group = "condition") + 
  expand_limits(y=0) +
  facet_wrap(~gene, nrow = 1, scales = "fixed", labeller = labeller(gene = qual.sensor.genes)) + 
  scale_x_discrete(labels = c("Control", "0.05 mM DLA", "0.5 mM DLA", "5 mM DLA", "45 mM DLA", "45 mM LLA")) + 
  theme(strip.text = element_text(face = "italic"), legend.position = 'none') + 
  scale_color_manual(values = c("black", blues9[5:8], "yellow4"))
p3d

p3 <- ggarrange(p3ab, p3c, p3d, nrow=3, labels = c("", "c", "d"), heights = c(2,3,2))
p3

png("figs/Fig3.png", width = 20, height = 25.3, units="cm", res=300)
p3
#ggsave("figs/Fig3.png", width = 20, height = 22, units = "cm") ##this should have worked, but fonts... :(
dev.off()


## just in case: all 10 DEGs for DLvsLL
myDegPlotWide(genes=c("YHL028W", "YLR054C", "YLR121C", "YGR189C", "YGR146C", 
                      "YGL255W", "YHR209W", "YLR205C", "YKR091W", "YCR005C"), counts = dds, group = "condition") + 
  #geom_boxplot() + 
  expand_limits(y=0) +
  facet_wrap(~gene, nrow = 2, scales = "fixed") + ##labeller = labeller(gene = gene_names)
  scale_x_discrete(labels = c("Control", "0.05 mM DLA", "0.5 mM DLA", "5 mM DLA", "45 mM DLA", "45 mM LLA")) + 
  theme(strip.text = element_text(face = "italic"), legend.position = 'none') + 
  scale_color_manual(values = c("black", blues9[5:8], "yellow4"))
#ggsave("DLvsLL_freescale.png", width = 24, height = 16, units="cm")  

