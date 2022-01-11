install.packages("bigmemory", repos="http://R-Forge.R-project.org")
install.packages(c("nortest","minpack.lm","Hmisc")

install.packages("iCAMP_1.3.1.tar.gz",repos = NULL, type = "source")
install.packages("NST_3.0.4.tar.gz",repos = NULL, type = "source")

library(iCAMP)

# version 2020.8.23
# version 2020.9.21 add classification information
rm(list=ls())
t0=Sys.time() # to calculate time cost

# 1 # set folder paths and file names, please change according to the folder paths and file names in your computer.
# the folder saving the input files
wd="/bigdata/wangmm/Stanford"

# the OTU table file (Tab delimited txt file)
com.file="OTU_cut1_no annotation.csv"

# the phylogenetic tree file
tree.file="Galaxy3-[FastTree.nwk].nwk"

# the classification (taxonomy) information
clas.file="OTU_Classifier_QIIME2_0.5.txt"

# the treatment informaiton table
treat.file="JRGCE_trt_water_addPlot1234.csv"

# the environmental varialbes
env.file="env_less.csv"
# if you do not have env file or the env may not represent niche, skip step 7 and 8, but check the alternative way to determine binning setting, e.g. bin.size.limit.

# the folder to save the output. please change to a new folder even if you are just testing the example data.
save.wd="/bigdata/wangmm/Stanford/iCAMP"

if(!dir.exists(save.wd)){dir.create(save.wd)}

# tNST
tnstout.jac=NST::tNST(comm=comm, group=treat.use, dist.method="jaccard", 
                  abundance.weighted = FALSE, rand=rand.time,  
                  nworker=nworker, null.model="PF", output.rand = TRUE, 
                  SES = FALSE, RC = FALSE)
write.csv(tnstout.jac$index.grp,file = paste0(prefix,".tNST.jac.summary.",colnames(treat)[i],".csv"))
write.csv(tnstout.jac$index.pair.grp,file = paste0(prefix,".tNST.jac.pairwise.",colnames(treat)[i],".csv"))


# 14.1b # bootstrapping test for tNST
tnstbt.boot.jac=nst.boot(nst.result=tnstout.jac, group=treat.use, rand=rand.time, trace=TRUE,
                two.tail=FALSE, out.detail=TRUE, between.group=FALSE, nworker=nworker)

save(tnstbt.boot.jac,file = paste0(prefix,".tNST.jac.boot.rda"))

write.table(tnst.boot.jac$summary,file = paste0(prefix,".tNST.jac.bootstr.",colnames(treat)[i],".csv"), quote = FALSE,sep = ",")
write.table(tnst.boot.jac$compare,file = paste0(prefix,".tNST.jac.compare.",colnames(treat)[i],".csv"), quote = FALSE,sep = ",")


##plot
library(usedist)
library(dplyr)
library(tibble)
library(reshape2)


##NST PERMANOVA
tNST <- read.csv(file = 'results//QPEN//STF.tNST.pairwise.water.csv', header = T, sep = ",")
dim(tNST)
setwd("D:\\lab_Tsinghua\\aProjects\\aStanford\\water_32 vs 32")

load(file = 'results\\QPEN\\tnstout.rda')
load(file = 'results\\QPEN\\pnstout.rda')
load(file = 'results\\QPEN\\tnstout.jac.rda')
head(tnstout$index.pair)
head(tnstout1$index.grp)
head(tnstout$index.pair.grp)
head(tnstout.jac$index.grp)
head(tnstout.jac$index.pair.grp)

head(pnstout$index.pair)
head(pnstout$index.grp)

tnst.all <- tnstout.jac$index.pair
dim(tnst.all)

#convert num into 'dist' num class
nams <- with (tnst.all, unique(c(as.character(name2), as.character(name1))))
ST.dist <- with(tnst.all, structure(ST.ij.jaccard,
                                      Size = length(nams),
                                      Labels = nams,
                                      Diag = FALSE,
                                      Upper = FALSE,
                                      Method = "user",
                                      class = "dist"))
length(nams)
Group<-t[match(nams,rownames(t)),]
dim(Group)

##ANOVA
adonis(ST.dist ~ water + plot, data=Group, permutations = 999)

MST.dist <- with(tnst.all, structure(MST.ij.jaccard,
                                    Size = length(nams),
                                    Labels = nams,
                                    Diag = FALSE,
                                    Upper = FALSE,
                                    Method = "user",
                                    class = "dist"))

adonis(MST.dist ~ water + plot, data=Group, permutations = 999)


##figure for tNST

##grouped bar plot tNST.jaccard and pNST.unw
#extract distance for aP and eP
##tNST.jac
tnst.jac.res <- dist_groups(MST.dist, Group$water)
head(tnst.jac.res)

dis.res.input <- tnst.jac.res %>%
  rename(tNST = Distance) %>%
  subset(Label == c('Within am', 'Within el')) %>%
  summarise(mean=mean(stochasticity.ratio), sd = sd(stochasticity.ratio)) %>%
  mutate(mean = mean*100, sd = sd*100)
  

plot.mst <- ggplot(dis.res.input, aes(x=Label, y=mean, fill= Label))+
  geom_bar(position = "dodge", stat = "identity", width = 0.5) +
  scale_fill_discrete(labels = c ("aP", "eP")) + #edit legend labels
  scale_x_discrete(labels = c ("aP", "eP"))+
  geom_errorbar(aes(ymin = mean-sd, ymax = mean + sd), colour="black", position = position_dodge(0.5), width=0.2, size=0.10) +
  labs(y="Estimated stochasticity (%)", x=NULL, fill="Treatment")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=9, colour = "black"),
        axis.text.y = element_text(size=9,colour = "black"))+
  geom_signif(y_position = 95, xmin=0.8, xmax=1.2, tip_length = 0,
              annotations = paste("italic(p) == 0.011"), parse = TRUE) + #add p value, p was parsed to plotmath
  ylim(0,100)

ggsave(plot.mst, filename = "results/QPEN/Stochasticity.ratio.pdf", width = 85, height = 60, units = "mm")


