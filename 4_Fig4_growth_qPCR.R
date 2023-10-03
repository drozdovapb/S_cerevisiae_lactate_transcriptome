library(openxlsx)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(ggplot2)
library(readxl)

## first, growth rate
dat <- read.xlsx("data/TableS2_LA_yeast_growth.xlsx", sheet="Data_pH_compensation") 

min(dat$OD600_init)
max(dat$OD600_init)

## how many repeats?
unique(dat$sample_id)
## within how many days?
unique(substr(dat$sample_id, 1, 8))


## different months to compare result
#dat$month <- substr(dat$sample_id, 1, 6)

dat <- dat[dat$condition %in% c("control", "DLA4500", "DLS4500"), ]



testres <- pairwise.wilcox.test(dat$RGR2, dat$condition, p.adjust.method = "holm", paired = TRUE)

## calculate statistics
stat.df <- data.frame(
  xmin = c(colnames(testres$p.value)[1], colnames(testres$p.value)[1],
            colnames(testres$p.value)[2]),
  xmax = c(rownames(testres$p.value)[1], rownames(testres$p.value)[2],
          rownames(testres$p.value)[2]),
  pvalue = c(testres$p.value[1,1], testres$p.value[2,1], testres$p.value[2,2]),
  annotation = c("*", "ns", "*"),
  y = c(.23, .2, .24))

pA <- 
ggplot(dat, aes(x=condition, y=RGR2, col = condition)) + 
  geom_boxplot(aes(fill=condition)) + 
  geom_jitter(width = 0.05) + 
  xlab("") + ylab("Relative growth rate") + 
  theme_bw() + 
  expand_limits(y=.25) + 
  scale_color_manual(values = c("black", blues9[7], blues9[7])) + 
  scale_fill_manual(values = c("grey", "pink", "grey")) + 
  theme(legend.position = 'none') + 
  scale_x_discrete(labels = c("Control", "45 mM \n DLA", "45 mM \n DLS")) + 
  geom_signif(inherit.aes = FALSE, data=stat.df,
                mapping = aes(annotations = annotation, xmin=xmin, xmax=xmax, y_position=y),
              manual=TRUE, map_signif_level = TRUE)
pA

### and qPCR
qpcr.dat <- read_xlsx("data/TableS5_by4742_all_repeats_qPCR_data.xlsx", sheet = "Results")

## check NTC and no-RT controls
table(qpcr.dat[qpcr.dat$Task == "NTC", "Cт"])
table(qpcr.dat[startsWith(qpcr.dat$`Sample Name`, "RT"), "Cт"])
## all undetermined or > 33. pass.

## remove control wells
qpcr.dat <- qpcr.dat[!(qpcr.dat$Task == "NTC"), ]
qpcr.dat <- qpcr.dat[!startsWith(qpcr.dat$`Sample Name`, "RT"), ]
## remove the group we don't need
qpcr.dat[!(grepl("DLA5", qpcr.dat$`Sample Name`)), ] -> qpcr.dat


## check Ct variance (sd)
table(qpcr.dat$`Cт SD` < 1)
#qpcr.dat[which(qpcr.dat$`Cт SD` > 1), ]
## remove non-amplified wells!
qpcr.dat <- qpcr.dat[qpcr.dat$Cт < 26, ]

## get reference genes
qpcr.dat <- aggregate(`Cт Mean` ~ `Sample Name` + `Target Name`, qpcr.dat, median)
names(qpcr.dat) <- c("Sample", "Target", "Cq")
library(tidyr)
qpcr.dat %>% pivot_wider(names_from = Target, values_from = Cq) -> qpcr.dat.wide
qpcr.dat.wide$refs <- sqrt(qpcr.dat.wide$ACT1*qpcr.dat.wide$CDC19)

qpcr.dat.wide %>% pivot_longer(cols=c("AQR1", "DLD3", "FIT2", "YPS3"), names_to="Target", values_to="Cq") -> qpcr.dat.wref
qpcr.dat.wref[!is.na(qpcr.dat.wref$Cq), ] -> qpcr.dat.wref
qpcr.dat.wref$deltaCq <- qpcr.dat.wref$refs- qpcr.dat.wref$Cq
qpcr.dat.wref$Group <- ""
qpcr.dat.wref$Group[grepl("K", qpcr.dat.wref$Sample)] <- "control"
qpcr.dat.wref$Group[grepl("DLS45", qpcr.dat.wref$Sample)] <- "DLS45"


gene_names <- c("YNL065W, AQR1", "YEL071W, DLD3", "YOR382W, FIT2", "YPL280W, YPS3")
names(gene_names) <- c("AQR1", "DLD3", "FIT2", "YPS3")


pB <- 
ggplot(qpcr.dat.wref, aes(Group, deltaCq, col = Group)) + 
  facet_wrap(~Target, nrow = 1, labeller = labeller(Target = gene_names), scales = "fixed") + 
  geom_boxplot() +  geom_point() + 
  scale_x_discrete(labels = c("Control", "45 mM \n DLS")) + 
  theme_bw() + 
  ylab("ΔCт") + xlab("") +
  scale_color_manual(values = c("black", blues9[7])) + 
  theme(strip.text = element_text(face = "italic"), legend.position = 'none',
        axis.text.x = element_text(angle = 0))


pB


ggarrange(pA, pB, widths = c(1, 2), labels = c("a", "b"))


png("figs/Fig4.png", width = 20, height = 8, units = "cm", res=300)
ggarrange(pA, pB, widths = c(1, 2), labels = c("a", "b"))
dev.off()
