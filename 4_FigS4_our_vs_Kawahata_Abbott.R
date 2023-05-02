library(ggplot2)
library(openxlsx)
library(ggpubr)

#### Kawahata, acid shock
## these data were taken from the supplement => to xlsx manually
kawahata <- read.xlsx(xlsxFile = "data/Kawahata_0.3percentLL_30min_DE.xlsx")
kawahata$Fold.change <- as.numeric(kawahata$Fold.change)
kawahata$logFC <- log2(kawahata$Fold.change)

## our to comparison
ours <- read.xlsx(xlsxFile = "data/DL_LL_DE_only.xlsx", sheet = "LL 45 mM")
names(ours)[1] <- "ORF"


## combine
combined <- merge(kawahata, ours, by = "ORF")

pS4A <- 
  ggplot(combined, aes(logFC, log2FoldChange)) + 
  xlab("Kawahata et al., 2006, acid shock") + ylab("Our data, 45 mM LLA") + 
  ylim(-3, 3) + xlim(-3, 3) + 
  geom_point() + #geom_smooth(method = "lm") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw()
#ggsave("ours_vs_Kawahata_DE_5genes.png")


## acid adaptation
kawahata.adapt <- read.xlsx(xlsxFile = "data/Kawahata_0.3percentLL_0.1to1_DE.xlsx")
kawahata.adapt$Fold.change <- as.numeric(kawahata.adapt$Fold.change)
kawahata.adapt$logFC <- log2(kawahata.adapt$Fold.change)

combined <- merge(kawahata.adapt, ours, by = "ORF")

pS4B <- 
ggplot(combined, aes(logFC, log2FoldChange)) + 
  xlab("Kawahata et al., 2006, acid adaptation") + ylab("Our data, 45 mM LLA") + 
  #  ylim(-3, 3) + xlim(-3, 3) + 
  geom_point() + #geom_smooth(method = "lm") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw()
pS4B
#ggsave("ours_vs_Kawahata.adapt_DE_5genes.png")



## now the same for Abbott, for which we had redone the calculations (see the previous script)


abbott500 <- read.delim("data/Abbott_DE_results_500.tsv")
## only DE
abbott500 <- abbott500[abbott500$adj.P.Val < 0.05 & abs(abbott500$logFC) > 1, ]
## our
ours <- read.xlsx(xlsxFile = "data/DL_LL_DE_only.xlsx", sheet = "LL 45 mM")
names(ours)[1] <- "ORF"
names(abbott500) <- gsub("Platform_ORF", "ORF", names(abbott500))


combined <- merge(abbott500, ours, by = "ORF")
## basically, there are just ours
pS4C <- 
ggplot(combined, aes(logFC, log2FoldChange)) + 
  ylab("Our data, 45 mM LLA") + xlab("Abbott et al., 2008, 500 mM LLA at pH 3") + 
  geom_point() + geom_smooth(method = "lm") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw()
pS4C
#ggsave("Abbot 500 mM DE only.png")

summary(lm(formula = log2FoldChange ~ logFC, data = combined))
##(Intercept) -0.38293    0.10546  -3.631 0.000946 ***
##  logFC        0.57087    0.04888  11.680 2.89e-13 ***
cor(combined$log2FoldChange, combined$logFC, use = "complete.obs")
## 0.89



abbott900 <- read.delim("data/Abbott_DE_results_900.tsv")
abbott900 <- abbott900[abbott900$adj.P.Val < 0.05 & abs(abbott900$logFC) > 1, ]
names(abbott900) <- gsub("Platform_ORF", "ORF", names(abbott900))
combined <- merge(abbott900, ours, by = "ORF")
## basically, there are just ours

pS4D <- 
  ggplot(combined, aes(logFC, log2FoldChange)) + 
  ylab("Our data, 45 mM LLA") + xlab("Abbott et al., 2008, 900 mM LLA at pH 5") + 
  geom_point() + geom_smooth(method = "lm") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw()
pS4D
#ggsave("Abbot 900 mM DE only.png")

summary(lm(formula = log2FoldChange ~ logFC, data = combined))
##(Intercept) -0.38293    0.10546  -3.631 0.000946 ***
##  logFC        0.57087    0.04888  11.680 2.89e-13 ***
cor(combined$log2FoldChange, combined$logFC, use = "complete.obs")
## 0.89, much worse

pS4 <- ggarrange(pS4A, pS4B, pS4C, pS4D, labels=LETTERS[1:4])
pS4

png("figs/FigS4.png", width=20, height=16, units="cm", res=300)
pS4
dev.off()
