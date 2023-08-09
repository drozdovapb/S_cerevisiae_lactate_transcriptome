library(openxlsx)
library(ggplot2)
library(ggbeeswarm)

dat <- read.xlsx("data/LA_yeast_growth.xlsx") ## rename S...

## doubling time doesn't work
#ggplot(dat, aes(x=condition, y=Doubling.time)) + 
#  geom_point() + ylim(0,10)

#a <- 
#ggplot(dat, aes(x=condition, y=Log2.Growth, col = condition)) + 
#  geom_boxplot() + 
#  geom_beeswarm() + 
#  #geom_jitter(width =.1) + 
#  xlab("") + ylab("log2 growth") + 
#  theme_bw() + 
#  scale_x_discrete(labels = c("control", "0.05 mM DLA", "0.5 mM DLA", "5 mM DLA", "45 mM DLA", "45 mM DLA")) + 
#  scale_color_manual(values = c("black", blues9[5:8], "yellow4"))
##ggsave("growth_inhibition.png", width = 20, height = 10, units="cm")

pairwise.wilcox.test(dat$Log2.Growth, dat$condition)


## relative growth rate!

png("figs/FigS1_growth_inhibition.png", width = 20, height = 8, units="cm", res=200)
ggplot(dat, aes(x=condition, y=RGR2, col = condition)) + 
  geom_boxplot() + 
  geom_beeswarm() + 
  #geom_jitter(width =.1) + 
  xlab("") + ylab("Relative growth rate") + 
  theme_bw() + 
  scale_x_discrete(labels = c("control", "0.05 mM DLA", "0.5 mM DLA", "5 mM DLA", "45 mM DLA", "45 mM LLA")) + 
  scale_color_manual(values = c("black", blues9[5:8], "yellow4")) + 
  theme(legend.position = 'none')
#ggsave("figs/FigS1_growth_inhibition.png", width = 20, height = 8, units="cm")
dev.off()

pairwise.wilcox.test(dat$RGR, dat$condition, p.adjust.method = "BH")
library(dunn.test)
dunn.test(x = dat$RGR, g = dat$condition, method = "Holm")

min(dat$OD600_init)
max(dat$OD600_init)
