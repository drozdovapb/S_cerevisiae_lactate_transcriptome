#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(openxlsx) ### data input
library(DESeq2) ## diff expression
library(EnhancedVolcano) ## plots Fig1
library(ggplot2) ## generally for plotting
library(ggpubr) ## for ggarrange


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

write.xlsx(x = list(res_DL0.05vsC.df, res_DL0.5vsC.df, res_DL5vsC.df, res_DL45vsC.df, res_LL45vsC.df, res_DL45vsLL45.df),
           sheetName = c("DL 0.05 mM", "DL 0.5 mM", "DL 5 mM", "DL 45 mM", "LL 45 mM", "DL vs LL"),
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
                      col = c("grey30", "grey30", "grey30", "red2"),
                      ylab="",
                      caption="", selectLab = "", legendPosition = 'none')
pD


pE <- EnhancedVolcano(res_DL0.05vsC, lab = rownames(res),
                      x = 'log2FoldChange', y = 'pvalue', 
                      pCutoff=0.05, pCutoffCol = 'padj', FCcutoff = 1,
                      title="e", subtitle="0.05 mM DLA vs. control",
                      col = c("grey30", "grey30", "grey30", "red2"),
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


p1 <- ggarrange(pA, pB, pC, pD, pE, pF)#, labels = letters[1:6])
png("figs/Fig1.png", width = 14, height = 10, units = "in", res=300)
#ggsave("figs/Fig1.png", width = 14, height = 10)
p1
dev.off()

save.image(file='data/DE.RData')

