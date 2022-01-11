if(!require(devtools))
  install.packages("devtools")
if(!require(ggcor))
  devtools::install_git("https://gitee.com/dr_yingli/ggcor")
library(vegan)
library(ape)
library(picante)
library(phyloseq)

library(ggplot2)
library(ggcor)
library(dplyr)

setwd("D:\\lab_Tsinghua\\aProjects\\aStanford\\water_32 vs 32")
ev <- read.table('env_less.csv', header=T, sep=",", row.names = 1)
dim(ev)
EV <- scale(ev)

#Mantel test
##For decreased OTUs
OTU.dec <- read.table('OTU_morethan0.01%_DESeq2&split-plot ANOVA_sig.dec.csv',header = T, sep = ",", row.names = 1)
dis.dec <- vegdist(t(OTU.dec))

m.dec <- matrix(0, nrow=24, ncol=2)
for (i in 1:24) {
  m.dec[i,1]<-mantel(dis.dec, dist(EV[,i],"euclid"), permutations = 999)$statistic;
  m.dec[i,2]<-mantel(dis.dec, dist(X2[,i],"euclid"), permutations = 999)$signif;
}
m.dec

##For increased OTUs
OTU.inc <- read.table('OTU_morethan0.01%_DESeq2&split-plot ANOVA_sig.inc.csv', header = T, sep = ",", row.names = 1)
dis.inc <- vegdist(t(OTU.inc))
m.inc <- matrix(0, nrow=24, ncol=2)
for (i in 1:24) {
  m.inc[i,1]<-mantel(dis.inc, dist(EV[,i],"euclid"), permutations = 999)$statistic;
  m.inc[i,2]<-mantel(dis.inc, dist(X2[,i],"euclid"), permutations = 999)$signif;
}
m.inc

##For stress genes
stress <- read.table('GeoChip_lnfill0mr_stress_sigchg.csv', header = T, sep = ",", row.names = 1)
dis.str <- vegdist(t(stress))
m.str <- matrix(0, nrow=24, ncol=2)
for (i in 1:24) {
  m.str[i,1]<-mantel(dis.str, dist(EV[,i],"euclid"), permutations = 999)$statistic;
  m.str[i,2]<-mantel(dis.str, dist(X2[,i],"euclid"), permutations = 999)$signif;
}
m.str

##Prepare mantel results for ggcor
mantel.result <- read.table('results\\mantel test_sigchg OTUs_for ggcor.csv', header = T, sep = ",")
dim(mantel.result)
head(mantel.result)
#convert to tibble

mantel.tibble <- as_tibble(mantel.result)
str(mantel.tibble)

mantel.sig <- mantel.tibble %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.1, 0.2, Inf),
                  labels = c("< 0.1", "0.1 - 0.2", "> 0.2")),
         pd = cut(p.value, breaks = c(-Inf, 0.05, 0.10, Inf),
                  labels = c("< 0.05", "0.05-0.10", "> 0.10")))

set_scale()

##change EV order
EV1 <- EV[,c(2:16,1)]
head(EV1)
p1.sig <- quickcor(EV1, type = "upper", grid.colour="grey",grid.size=0.5, show.diag = FALSE) +
  geom_colour(width=0.7, height=0.7, colour = "white") +
  anno_link(data = mantel.sig, aes(colour=pd, size=rd)) +
  scale_size_manual(values = c(0.75, 1, 2)) +
  scale_colour_manual(values = c("#993399", "#1B9E77", "#D2D2D2")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), reverse = TRUE,
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 2), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))

ggsave(p1.sig, filename = "results//ggcor.sig.pdf")
