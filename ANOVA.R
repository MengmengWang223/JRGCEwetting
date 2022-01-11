setwd("D:\\lab_Tsinghua\\aProjects\\aStanford\\water_32 vs 32")

##split-plot ANOVA

library(vegan)
library(agricolae)
library(ggprism)
#for environmental variables, with samples in rows and variables in columns
All<-read.table('env_less.csv', header=T, sep=",", row.names = 1)
#for diversity
All<-read.table('results\\alpha diversity.csv', header=T, sep=",", row.names = 1)
#for GeoChip category and subcategory
All<-read.table('GeoChip_lnfill0mr_category and specific subcategory.csv', header=T, sep=",", row.names = 1)
All <- t(All)
##for GeoChip gene
All<-read.table('GeoChip_lnfill0mr_summed to gene_phage.csv', header=T, sep=",", row.names = 1)
All<-read.table('GeoChip_lnfill0mr_summed to gene_nitrogen.csv', header=T, sep=",", row.names = 1)
All<-read.table('GeoChip_lnfill0mr_summed_virus subcategory&gene.csv', header=T, sep=",")
All <- t(All)
##for phylum summary
All <- read.table('OTU_phylum.comb summary.csv', header=T, sep=",", row.names = 1)
All <- read.table('OTU_class.comb summary.csv', header=T, sep=",")
All <- read.table('OTU_comb summary_genus above_morethan10.csv', header=T, sep=",")
dim(All)
#All <- t(All[,-c(1:5)])

All <- t(All[,-1])
head(All)
dim(All)
attributes(All) 

Y<-data.matrix(All, rownames.force = NA)
head(Y)


Group=read.table('D:\\lab_Tsinghua\\aProjects\\aStanford\\water_32 vs 32\\JRGCE_trt_water_addPlot1234.csv',header=T, sep=",", row.names = 1)
head(Group)
Group1<-Group[match(rownames(All),rownames(Group)),]
head(Group1)
dim(Group1)


#randomized block split-plot anova
aov.ep=aov(Y ~ water+Error(block/plot1234), data=Group1)
summary(aov.ep)

aov.ex2=aov(Y ~ CO2*heat*water*nitrogen+Error(block/plot1234), data=Group1)
summary(aov.ex2)
tab<-summary(aov.ex2)

##export the aov results to a text file rather than the console
sink("results\\ANOVA_GeoChip_ln0mr_virus subcategory&genes.txt")
print(summary(aov.ep))
sink() # returns output to the console


###split-plot for each OTU. Since the OTU number is too large, we just extract main water effect on each OTU rather than 
##the whole results of split-plot ANOVA.
library(broom)
OTU <- read.table("OTU_nosinglet_morethan0.05%.csv", header = TRUE, sep = ",", row.names = 1)
dim(OTU)
OTU.anno <- OTU[,1:10]
OTU.data <- OTU[,-(1:10)]
dim(OTU.data)

Group=read.table('D:\\lab_Tsinghua\\aProjects\\aStanford\\water_32 vs 32\\JRGCE_trt_water_addPlot1234.csv',header=T, sep=",", row.names = 1)
head(Group)
Group1<-Group[match(colnames(OTU.data),rownames(Group)),]
head(Group1)

## do ANOVA for each OTU and extract p value of water effect
library(rstatix)
length = nrow(OTU.data)
p <- matrix(0, nrow=length, ncol=1)
  
for(i in 1:length){
    Y <- t(data.matrix(OTU.data[i,], rownames.force = NA))
    aov.sp=aov(Y ~ water+Error(block/plot1234), data=Group1)
    tab <- anova_summary(aov.sp)
    p[i] <- tab$"p"[1]
    }
p.adjusted <- apply(p,2,function(x)(p.adjust(x, method = "BH", n = length(x))))
p.adjust(p,method="BH", n=length(p))
res.p <- bind_cols(OTU.anno, p.value = p, p.adjusted = p.adjusted)
write.csv(res.p, file = "results/OTU&DESeq2&split-plot ANOVA/OTU_morethan0.05%_split-plot ANOVA results.csv")
