library(openxlsx)
library(ggplot2)
library(ggbeeswarm)

## read data
dat <- read.xlsx("data/TableS2_LA_yeast_growth.xlsx", sheet="Data_growth_inhibition") 

## check if there are statistically significant changes
library(dunn.test)
dunn.test(x = dat$RGR2, g = dat$condition, method = "Holm")

## check min and max values for initial OD600
min(dat$OD600_init)
max(dat$OD600_init)

## how many repeats?
unique(dat$sample_id)
## within how many days?
unique(substr(dat$sample_id, 1, 8))

pS1 <- 
ggplot(dat, aes(x=condition, y=RGR2, col = condition)) + 
  geom_boxplot() + 
  geom_beeswarm() + 
  xlab("") + ylab("Relative growth rate") + 
  theme_bw() + 
  scale_x_discrete(labels = c("control", "0.05 mM DLA", "0.5 mM DLA", 
                              "5 mM DLA", "45 mM DLA", "45 mM LLA")) + 
  scale_color_manual(values = c("black", blues9[5:8], "yellow4")) + 
  theme(legend.position = 'none')

png("figs/FigS1_growth_inhibition.png", width = 20, height = 8, units="cm", res=200)
pS1
dev.off()
