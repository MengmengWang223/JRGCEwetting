library(vegan)

#import data
Seq=read.table('D:\\lab_Tsinghua\\aProjects\\aStanford\\water_32 vs 32\\OTU_nosinglet_no annotation.csv',header=T,sep=",", row.names = 1)
dim(Seq)
# or
Seq=read.table('OTU_no_singlet_no_mitochondria_no_chloroplast_noanno_resample26658.csv',header=T,sep=",", row.names = 1)
dim(Seq)
#
Geo=read.table('D:\\lab_Tsinghua\\aProjects\\aStanford\\water_32 vs 32\\GeoChip_fill0lnmr_no annotation.csv',header=T,sep=",")
dim(Geo)

t=read.table('D:\\lab_Tsinghua\\aProjects\\aStanford\\water_32 vs 32\\JRGCE_trt_water_addPlot1234.csv',header=T,sep=",", row.names = 1)

head(t)
Group<-t[match(colnames(Seq),rownames(t)),]
head(Group)
dim(Group)
str(Group)

##4 aP with 4 eP
Seq.p = Seq[,c(1:4, 33:36)]
dim(Seq.p)
Group.p = Group[c(1:4, 33:36),]
Group.p
seq.p.dist <- vegdist(t(Seq.p), method = "bray")

adonis(seq.p.dist ~ water, data=Group.p, permutations=999)

##ANOVA with split-plot error
adonis(t(Seq) ~ water + plot, data=Group, permutations = 999)

##split-plot PERMANOVA
library(usedist)
library(dplyr)
library(tibble)
seq.dist <- vegdist(t(Seq), method = "bray")

###Top-half: on plot level, for block random error and plot-level factor
wholeplot <- factor(Group$plot)

# to get the full centroid matrix
seq.centroid.1 <- dist_multi_centroids(seq.dist, wholeplot, squared = F)
seq.centroid.matrix<-seq.centroid$vectors
dim(seq.centroid.matrix)

## rearrange CO2 and heat treatment on plot level
Group.plot <- Group %>% 
  select(plot, CO2, heat, block) %>%
  distinct(plot, .keep_all = T) %>%
  arrange(plot) %>% remove_rownames %>% 
  column_to_rownames(var = "plot")

dim(Group.plot)

adonis(seq.centroid.1 ~ block + CO2*heat, data=Group.plot, strata=Group.plot$block, method = "euclidean", permutations = 9999)
adonis(seq.centroid.1 ~ block + CO2*heat, data=Group.plot, method = "euclidean", permutations = 9999)

#Lower-half: for sub-plot factor,whole-plot*sub-plot interaction, plot error
adonis(t(Seq) ~ water*nitrogen + CO2:water + CO2:nitrogen + CO2:water:nitrogen + heat:water + heat:nitrogen 
       + heat:water:nitrogen + CO2:heat:water + CO2:heat:nitrogen + heat:CO2:water:nitrogen + plot, data=Group, permutations = 999)
