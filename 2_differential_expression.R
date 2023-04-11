#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)

count_table <- read.delim("allSamples.featureCounts.txt", skip=1, row.names="Geneid")
sample_table <- read.delim("sample_table.txt", row.names = "sample")


dds <- DESeqDataSetFromMatrix(countData = count_table[,6:23], colData = sample_table, design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
#write.table(as.data.frame(resOrdered),sep='\t', quote=FALSE, file='wt_mutant_p-values.txt')

vsd <- vst(dds)
plotPCA(vsd, intgroup=c("condition"))


res_DL45vsC = results(dds, contrast=c("condition", "DL4500", "control"))
plot <- plotMA(res_DL45vsC, main = 'DL45', ylim = c(-2,2), xlab = 'mean count')
res_DL45vsC.df <- as.data.frame(res_DL45vsC)
res_DL45vsC.df <- res_DL45vsC.df[order(res_DL45vsC.df$padj), ]

res_DL5vsC = results(dds, contrast=c("condition", "DL0500","control"))
plot <- plotMA(res_DL5vsC, main = 'DL5', ylim = c(-2,2), xlab = 'mean count')
res_DL5vsC.df <- as.data.frame(res_DL5vsC)
res_DL5vsC.df <- res_DL5vsC.df[order(res_DL5vsC.df$padj), ]

res_DL0.5vsC = results(dds, contrast=c("condition", "DL0050","control"))
plot <- plotMA(res_DL0.5vsC, main = 'DL0.5', ylim = c(-2,2), xlab = 'mean count')
res_DL0.5vsC.df <- as.data.frame(res_DL0.5vsC)
res_DL0.5vsC.df <- res_DL0.5vsC.df[order(res_DL0.5vsC.df$padj), ]

res_DL0.05vsC = results(dds, contrast=c("condition", "DL0005","control"))
plot <- plotMA(res_DL0.05vsC, main = 'DL0.05', ylim = c(-2,2), xlab = 'mean count')
res_DL0.05vsC.df <- as.data.frame(res_DL0.05vsC)
res_DL0.05vsC.df <- res_DL0.05vsC.df[order(res_DL0.05vsC.df$padj), ]

res_LL45vsC = results(dds, contrast=c("condition", "LL4500", "control"))
plot <- plotMA(res_LL45vsC, main = 'LL45', ylim = c(-2,2), xlab = 'mean count')
res_LL45vsC.df <- as.data.frame(res_LL45vsC)
res_LL45vsC.df <- res_LL45vsC.df[order(res_LL45vsC.df$padj), ]

library(openxlsx)
write.xlsx(x = list(res_DL0.05vsC.df, res_DL0.5vsC.df, res_DL5vsC.df, res_DL45vsC.df, res_LL45vsC.df),
           sheetName = c("DL 50 uM", "DL 500 uM", "DL 5 mM", "DL 45 mM", "LL 45 mM"),
           asTable = T, rowNames = T, file = "DL_LL_DE.xlsx")



## and also differentially expressed genes
fcthreshold <- sqrt(1.5)
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

write.xlsx(x = list(DL0005.degs, DL0050.degs, DL0500.degs, DL4500.degs, LL4500.degs),
           sheetName = c("DL 50 uM", "DL 500 uM", "DL 5 mM", "DL 45 mM", "LL 45 mM"),
           asTable = T, rowNames = T, file = "DL_LL_DE_only.xlsx")


## heatmap for 133 top genes
topGenes <- head(order(res$padj), 10)
betas <- coef(dds)
mat <- betas[topGenes, -c(1)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
library("pheatmap")
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=TRUE)
