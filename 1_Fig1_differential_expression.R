#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(openxlsx) ### data input
library(DESeq2) ## diff expression
library(EnhancedVolcano) ## plots Fig1
library(ggplot2) ## generally for plotting
library(ggpubr) ## for ggarrange
library(BioVenn) ## for venn diagram with proportional circles
library(cowplot) ## for venn diagram with proportional circles

## read data
count_table <- read.delim("data/allSamples.featureCounts.txt", skip=1, row.names="Geneid")
sample_table <- read.delim("data/sample_table.txt", row.names = "sample")


dds <- DESeqDataSetFromMatrix(countData = count_table[,6:23], colData = sample_table, design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
#write.table(as.data.frame(resOrdered),sep='\t', quote=FALSE, file='wt_mutant_p-values.txt')

## save normalized counts for GEO submission
normalized_counts_DESeq <- counts(dds, normalized=TRUE)
write.csv(normalized_counts_DESeq, file="data/normalized_counts_DESeq.csv")

## check PCA (nothing thrilling)
vsd <- vst(dds)
plotPCA(vsd, intgroup=c("condition"))

## get results
res_DL45vsC = results(dds, contrast=c("condition", "DL4500", "control"))
res_DL45vsC.df <- as.data.frame(res_DL45vsC)
res_DL45vsC.df <- res_DL45vsC.df[order(res_DL45vsC.df$padj), ]

res_DL5vsC = results(dds, contrast=c("condition", "DL0500","control"))
res_DL5vsC.df <- as.data.frame(res_DL5vsC)
res_DL5vsC.df <- res_DL5vsC.df[order(res_DL5vsC.df$padj), ]

res_DL0.5vsC = results(dds, contrast=c("condition", "DL0050","control"))
res_DL0.5vsC.df <- as.data.frame(res_DL0.5vsC)
res_DL0.5vsC.df <- res_DL0.5vsC.df[order(res_DL0.5vsC.df$padj), ]

res_DL0.05vsC = results(dds, contrast=c("condition", "DL0005","control"))
res_DL0.05vsC.df <- as.data.frame(res_DL0.05vsC)
res_DL0.05vsC.df <- res_DL0.05vsC.df[order(res_DL0.05vsC.df$padj), ]

res_LL45vsC = results(dds, contrast=c("condition", "LL4500", "control"))
res_LL45vsC.df <- as.data.frame(res_LL45vsC)
res_LL45vsC.df <- res_LL45vsC.df[order(res_LL45vsC.df$padj), ]

res_DL45vsLL45 = results(dds, contrast=c("condition", "DL4500", "LL4500"))
res_DL45vsLL45.df <- as.data.frame(res_DL45vsLL45)
res_DL45vsLL45.df <- res_DL45vsLL45.df[order(res_DL45vsLL45.df$padj), ]

#### different concentrations 
## 45 mM vs 5 mM
res_DL45vsDL05 = results(dds, contrast=c("condition", "DL4500", "DL0500"))
res_DL45vsDL05.df <- as.data.frame(res_DL45vsDL05)
res_DL45vsDL05.df <- res_DL45vsDL05.df[order(res_DL45vsDL05.df$padj), ]
DL4505.degs <- 
res_DL45vsDL05.df[abs(res_DL45vsDL05.df$log2FoldChange) > fcthreshold & res_DL45vsDL05.df$padj < 0.05 & 
                    complete.cases(res_DL45vsDL05.df$padj), ]


#res_DL45vsDL0050 = results(dds, contrast=c("condition", "DL4500", "DL0050"))
#res_DL45vsDL0050.df <- as.data.frame(res_DL45vsDL0050)
#res_DL45vsDL0050.df <- res_DL45vsDL0050.df[order(res_DL45vsDL0050.df$padj), ]
#DL45_0050.degs <- 
#res_DL45vsDL0050.df[abs(res_DL45vsDL0050.df$log2FoldChange) > fcthreshold & res_DL45vsDL0050.df$padj < 0.05 & 
#                    complete.cases(res_DL45vsDL0050.df$padj), ]


#res_DL45vsDL0005 = results(dds, contrast=c("condition", "DL4500", "DL0005"))
#res_DL45vsDL0005.df <- as.data.frame(res_DL45vsDL0050)
#res_DL45vsDL0005.df <- res_DL45vsDL0005.df[order(res_DL45vsDL0005.df$padj), ]
#DL45_0005.degs <- 
#  res_DL45vsDL0005.df[abs(res_DL45vsDL0005.df$log2FoldChange) > fcthreshold & res_DL45vsDL0005.df$padj < 0.05 & 
#                        complete.cases(res_DL45vsDL0005.df$padj), ]

####

write.xlsx(x = list(res_DL0.05vsC.df, res_DL0.5vsC.df, res_DL5vsC.df, res_DL45vsC.df, res_LL45vsC.df, res_DL45vsLL45.df, res_DL45vsDL05.df),
           sheetName = c("DL 0.05 mM", "DL 0.5 mM", "DL 5 mM", "DL 45 mM", "LL 45 mM", "DL vs LL", "DL 45 vs 5 mM"),
           asTable = T, rowNames = T, file = "data/TableS1_DL_LL_DE.xlsx")



## and also differentially expressed genes
fcthreshold <- 1
DL0005.degs <- res_DL0.05vsC.df[abs(res_DL0.05vsC.df$log2FoldChange) > fcthreshold & res_DL0.05vsC.df$padj < 0.05 & 
                                  complete.cases(res_DL0.05vsC.df$padj), ]
DL0050.degs <- res_DL0.5vsC.df[abs(res_DL0.5vsC.df$log2FoldChange) > fcthreshold & res_DL0.5vsC.df$padj < 0.05 & 
                                 complete.cases(res_DL0.5vsC.df$padj), ]
DL0500.degs <- res_DL5vsC.df[abs(res_DL5vsC.df$log2FoldChange) > fcthreshold & res_DL5vsC.df$padj < 0.05 & 
                               complete.cases(res_DL5vsC.df$padj), ]
DL4500.degs <- res_DL45vsC.df[abs(res_DL45vsC.df$log2FoldChange) > fcthreshold & res_DL45vsC.df$padj < 0.05 & 
                                complete.cases(res_DL45vsC.df$padj), ]
LL4500.degs <- res_LL45vsC.df[abs(res_LL45vsC.df$log2FoldChange) > fcthreshold & res_LL45vsC.df$padj < 0.05 & 
                                complete.cases(res_LL45vsC.df$padj), ]
DLLL.degs <- res_DL45vsLL45.df[abs(res_DL45vsLL45.df$log2FoldChange) > fcthreshold & res_DL45vsLL45.df$padj < 0.05 & 
                                 complete.cases(res_DL45vsLL45.df$padj), ]
DL4505.degs <- res_DL45vsDL05.df[abs(res_DL45vsDL05.df$log2FoldChange) > fcthreshold & res_DL45vsDL05.df$padj < 0.05 & 
                                 complete.cases(res_DL45vsDL05.df$padj), ]

write.xlsx(x = list(DL0005.degs, DL0050.degs, DL0500.degs, DL4500.degs, LL4500.degs, DLLL.degs),
           sheetName = c("DL 50 uM", "DL 500 uM", "DL 5 mM", "DL 45 mM", "LL 45 mM", "DLvsLL"),
           asTable = T, rowNames = T, file = "data/DL_LL_DE_only.xlsx")


## heatmap for 133 top genes (not used in the final manuscript)
topGenes <- head(order(res$padj), 10)
betas <- coef(dds)
mat <- betas[topGenes, -c(1)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
library("pheatmap")
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=TRUE)


## make some volcanoes!
pA <- EnhancedVolcano(res_DL45vsC, lab = rownames(res),
                      x = 'log2FoldChange', y = 'pvalue', 
                      pCutoff=0.05, pCutoffCol = 'padj', FCcutoff = 1,
                      title="a", subtitle="45 mM DLA vs. control", 
                      col = c("grey30", "grey30", "grey30", "red2"),
                      xlab="",
                      ylab = bquote(~-Log[10] ~ italic(p)),
                      caption="", selectLab = "", legendPosition = 'none')
pA
pB <- EnhancedVolcano(res_DL5vsC, lab = rownames(res),
                      x = 'log2FoldChange', y = 'pvalue', 
                      pCutoff=0.05, pCutoffCol = 'padj', FCcutoff = 1,
                      title="b", subtitle="5 mM DLA vs. control",
                      col = c("grey30", "grey30", "grey30", "red2"),
                      xlab = "", ylab="",
                      caption="", selectLab = "", legendPosition = 'none')
pB

pC <- EnhancedVolcano(res_LL45vsC, lab = rownames(res),
                      x = 'log2FoldChange', y = 'pvalue', 
                      pCutoff=0.05, pCutoffCol = 'padj', FCcutoff = 1,
                      title="c", subtitle="45 mM LLA vs. control",
                      col = c("grey30", "grey30", "grey30", "red2"),
                      xlab="", ylab="",
                      caption="", selectLab = "", legendPosition = 'none')
pC


pD <- EnhancedVolcano(res_DL0.5vsC, lab = rownames(res),
                      x = 'log2FoldChange', y = 'pvalue', 
                      pCutoff=0.05, pCutoffCol = 'padj', FCcutoff = 1,
                      title="d", subtitle="0.5 mM DLA vs. control",
                      ylab = bquote(~-Log[10] ~ italic(p)),
                      col = c("grey30", "grey30", "grey30", "red2"),
                      caption="", selectLab = "", legendPosition = 'none')
pD


pE <- EnhancedVolcano(res_DL0.05vsC, lab = rownames(res),
                      x = 'log2FoldChange', y = 'pvalue', 
                      pCutoff=0.05, pCutoffCol = 'padj', FCcutoff = 1,
                      title="e", subtitle="0.05 mM DLA vs. control",
                      col = c("grey30", "grey30", "grey30", "red2"),
                      ylab="",
                      caption="", selectLab = "", legendPosition = 'none')
pE


pF <- EnhancedVolcano(res_DL45vsLL45, lab = rownames(res),
                      x = 'log2FoldChange', y = 'pvalue', 
                      title="f", subtitle = "45 mM DLA vs. 45 mM LLA",
                      pCutoffCol = 'padj', pCutoff = 0.05, 
                      ylab="",
                      col = c("grey30", "grey30", "grey30", "red2"),
                      selectLab = "", legendPosition = 'none')
pF


pNewA <- EnhancedVolcano(res_DL45vsLL45, lab = rownames(res),
                         x = 'log2FoldChange', y = 'pvalue', 
                         title="f", subtitle = "45 mM DLA vs. 45 mM LLA",
                         pCutoffCol = 'padj', pCutoff = 0.05, 
                         ylab="",
                         col = c("grey30", "grey30", "grey30", "red2"),
                         selectLab = "", legendPosition = 'none')

## venn diagram to see the intersections between the DEGs
DL5 = row.names(DL0500.degs) 
DL45 = row.names(DL4500.degs)
LL45 = row.names(LL4500.degs)
## upregs & downregs
DL5up <- row.names(DL0500.degs[DL0500.degs$log2FoldChange > 1, ])
DL5down <- row.names(DL0500.degs[DL0500.degs$log2FoldChange < -1, ])
DL45up <- row.names(DL4500.degs[DL4500.degs$log2FoldChange > 1, ])
DL45down <- row.names(DL4500.degs[DL4500.degs$log2FoldChange < -1, ])
LL45up <- row.names(LL4500.degs[LL4500.degs$log2FoldChange > 1, ])
LL45down <- row.names(LL4500.degs[LL4500.degs$log2FoldChange < -1, ])

save.image(file='data/DE.RData')

par(mar=c(0,2.5,4,2.5)) ## bottom, left, top, right default c(5, 4, 4, 2) + 0.1
draw.venn(DL5up, DL45up, LL45up, title = paste0("g", strrep(" ", 45), "\n \n"), 
          subtitle = "Upregulated genes", 
          xtitle = "5 mM DLA", ytitle = "45 mM DLA", ztitle = "45 mM LLA", 
          x_c = blues9[7], y_c=blues9[8], z_c="yellow4", 
          xt_f = "Arial", yt_f = "Arial", zt_f = "Arial", 
          xt_fb = 1, yt_fb = 1, zt_fb = 1, 
          xt_s = 1.2, yt_s = 1.2, zt_s = 1.2, nr_s = 1.2,
          nr_f = "Arial", t_f = "Arial", t_fb = 2, bg_c = "transparent",
          st_fb = 1, st_f = "Arial")
text(0.42, 0.975, "Upregulated genes", cex=1.2)
pG <- recordPlot()

par(mar=c(0,2.5,4,2.5)) ## bottom, left, top, right default c(5, 4, 4, 2) + 0.1
draw.venn(DL5down, DL45down, LL45down, title = paste0("h", strrep(" ", 45), "\n \n"), 
          subtitle="",
          xtitle = "5 mM DLA \n", ytitle = "45 mM DLA", ztitle = "\n 45 mM LLA", 
          x_c = blues9[7], y_c=blues9[8], z_c="yellow4", 
          xt_f = "Arial", yt_f = "Arial", zt_f = "Arial", 
          xt_fb = 1, yt_fb = 1, zt_fb = 1, 
          xt_s = 1.2, yt_s = 1.2, zt_s = 1.2, nr_s = 1.2,
          nr_f = "Arial", t_f = "Arial", t_fb = 2, bg_c = "transparent", 
          st_fb = 1, st_f = "Arial")
text(0.4, 0.975, "Downregulated genes", cex=1.2)
pH <- recordPlot()



#p1 <- ggarrange(pA, pB, pC, pD, pE, pF)#, labels = letters[1:6])
p1 <- plot_grid(plotlist = list(pA, pB, pC, pG, pD, pE, pF, pH), ncol = 4)
#p1

#library(cowplot)
#p1upd <- plot_grid(p1, pG, ncol=2, rel_widths = c(3, 1), greedy = FALSE)


#png("figs/Fig1.png", width = 14, height = 10, units = "in", res=300)
png("figs/Fig1.png", width = 18, height = 10, units = "in", res=300)
#ggsave("figs/Fig1.png", width = 14, height = 10)
p1
dev.off()


