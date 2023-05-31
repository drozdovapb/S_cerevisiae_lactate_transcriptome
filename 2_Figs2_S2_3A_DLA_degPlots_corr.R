load("DE.RData")

#library(ggplot2)
#library(ggpubr)

## Fig. 2: lactate-related genes
source("0_helper_functions.R")

gene_names <- c("YDL174C, DLD1", "YEL071W, DLD3", "YMR322C, SNO4", "YPL280W, HSP32")
names(gene_names) <- c("YDL174C", "YEL071W", "YMR322C", "YPL280W")

p2A <- 
  myDegPlotWide(dds, c("YDL174C", "YEL071W", "YMR322C", "YPL280W"), group="condition") + 
  expand_limits(y=c(0, 16)) + 
  facet_wrap(~gene, nrow = 1, labeller = labeller(gene = gene_names), scales = "fixed") + 
  scale_x_discrete(labels = c("control", "0.05 mM DLA", "0.5 mM DLA", "5 mM DLA", "45 mM DLA", "45 mM LLA")) + 
  theme(strip.text = element_text(face = "italic"), legend.position = 'none') + 
  scale_color_manual(values = c("black", blues9[5:8], "yellow4"))
p2A
#ggsave("lactate_metab.png", width = 20, height = 8, units="cm")  


## Fig. 2: correlation between DL45 and DL5
## different concentrations of DLA
#intersect(row.names(DL0500.degs), row.names(DL4500.degs))
DLmerged <- merge(DL0500.degs, DL4500.degs, by = 0)
cor(DLmerged$log2FoldChange.x, DLmerged$log2FoldChange.y)
pDLmerged <- 
  ggplot(DLmerged, aes(log2FoldChange.x, log2FoldChange.y)) + 
  geom_point(alpha=1) + theme_bw() + 
  xlab("log2 fold change, 5 mM DLA") + ylab("log2 fold change, 45 mM DLA") +
  geom_vline(xintercept=0, linetype="dotted") + geom_hline(yintercept=0, linetype="dotted")
pDLmerged

DLmergedAll <- merge(res_DL5vsC.df, res_DL45vsC.df, by = 0, all = TRUE)
cor(DLmergedAll$log2FoldChange.x, DLmergedAll$log2FoldChange.y, use="complete.obs")
pDLmergedAll <- 
  ggplot(DLmergedAll, aes(log2FoldChange.x, log2FoldChange.y)) + 
  geom_point(alpha=.1) + theme_bw() + 
  xlab("log2 fold change, 5 mM DLA") + ylab("log2 fold change, 45 mM DLA") +
  geom_vline(xintercept=0, linetype="dotted") + geom_hline(yintercept=0, linetype="dotted")
pDLmergedAll

p2 <- ggarrange(p2A, labels = "a", nrow=2,
                ggarrange(pDLmerged, pDLmergedAll, labels = c("b", "c")))
p2
#ggsave("figs/Fig2.png", width = 20, height = 15, units = "cm") ## should have worked but fonts... ;(
png("figs/Fig2.png", width = 20, height = 15, units = "cm", res=300)
p2
dev.off()
## and this is S2: cell wall (biosynthesis) genes

pS2 <- myDegPlotWide(genes=c("YAL053W", "YGR189C", "YMR305C", "YGR032W"), counts = dds, group = "condition") + 
  expand_limits(y=0) +
  #  facet_wrap(~gene, nrow = 1, scales = "fixed") + ##labeller = labeller(gene = gene_names)
  scale_x_discrete(labels = c("control", "50 uM DLA", "500 uM DLA", "5 mM DLA", "45 mM DLA", "45 mM LLA")) + 
  theme(strip.text = element_text(face = "italic"), legend.position = 'none') + 
  scale_color_manual(values = c("black", blues9[5:8], "yellow4"))
pS2

#ggsave("figs/Figs2_cell_wall.png", width = 20, height = 8, units="cm") ## should have worked but fonts... :(
png("figs/Figs2_cell_wall.png", width = 20, height = 8, units="cm", res=300)  
pS2
dev.off()

