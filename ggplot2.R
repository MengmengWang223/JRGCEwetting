setwd("D:\\lab_Tsinghua\\aProjects\\aStanford\\water_32 vs 32")
library(ggplot2)
library(ggprism)
library(patchwork)
library(magrittr)

data<-read.table(file = 'results//',header=T,sep=",")
dim(data)
head(data)

colnames(data)

data$Gene.name<-factor(data$Gene.name,levels=unique(data$Gene.name))
data$subcategory<-factor(data$subcategory,levels=unique(data$subcategory))
#stress gene bar plot with errorbar
bar<-ggplot(data=data,aes(x=Gene.name,y=Percentage))+
  geom_bar(stat = "identity",position="identity",fill="grey",alpha=0.7,color="black",size=0.05)+
  geom_errorbar(aes(ymin=Percentage-se,ymax=Percentage+se),colour="black",width=0.3,size=0.10)+
  ylim(0,22)+
  #theme(axis.text.x = element_text(angle))
  labs(x=NULL,y="Percent change (%)", size=10)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,vjust=1,size=9, colour = "black"),
        axis.text.y = element_text(size=9,colour = "black"))
bar
pdf("D:\\lab_Tsinghua\\aProjects\\aStanford\\water_32 vs 32\\results\\stress gene bar plot_2.pdf", width=5,height=3)
bar
dev.off()

##C degradation bar plot
data$Gene.name<-factor(data$Gene.name,levels=unique(data$Gene.name))
data$subcategory<-factor(data$subcategory,levels=unique(data$subcategory))


Var1=c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd")
names(Var1)=c("glucose","inulin","starch","pectin","agar","hemicellulose","cellulose","chitin","cutin","protein")

bar<-ggplot(data=data,aes(x=Gene.name,y=Percentage,fill=subcategory))+
       geom_bar(stat = "identity",position="identity",color="black",size=0.05)+
       geom_errorbar(aes(ymin=Percentage-se,ymax=Percentage+se),colour="black",width=0.25,size=0.10)+
       scale_fill_manual(values=Var1)+
    #theme(axis.text.x = element_text(angle))
    labs(x=NULL,y="Percent change (%)")+
    scale_y_continuous(breaks=seq(-5,20,5), limits=c(-9,22))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,vjust=1,size=9, colour = "black"),
        axis.text.y = element_text(size=9,colour = "black"))

pdf("D:\\lab_Tsinghua\\aProjects\\aStanford\\water_32 vs 32\\results\\Geochip\\fill0MR_Cdegradation bar_with errorbar and protein genes2.pdf", width=4.5,height=3)
bar
dev.off()

##barplot of environmental varaibles: soil moisture, total C, aboveground biomass& litter, shallow root & fine root

ev <- read.table('env_less.csv', header=T, sep=",", row.names = 1)
group <- read.table('JRGCE_trt_water.csv', header=T, sep=",", row.names = 1)

data.cb <- merge(ev, group, by="row.names", all.x=TRUE, sort = FALSE)
head(data.cb)

library(plotrix) ##containing se calculation
library(reshape2)
library(data.table)
library(dplyr)
ev.mean <- data.cb %>% group_by(water) %>% summarise_at(vars(Temperature:NH4), mean) %>% 
  reshape2::melt(id.vars = "water", variable.name = "ev", value.name = "mean")

           
ev.se <- data.cb %>% group_by(water) %>% summarise_at(vars(Temperature:NH4), std.error) %>%
  reshape2::melt(id.vars = "water", variable.name = "ev", value.name = "se")

data.input <- bind_cols(ev.mean, se=ev.se$se)

moisture <- data.input[data.input$ev == "Moisture",]

bar.sm <- ggplot(data=moisture,aes(x=ev, y=mean, fill=water))+
  geom_col(width =0.5, position = "dodge", size=0.05)+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),colour="black",position = position_dodge(0.5), width=0.2,size=0.10)+
  #theme(axis.text.x = element_text(angle))
  labs(x=NULL,y="Soil moisture (%)", size=10)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=9, colour = "black"),
        axis.text.y = element_text(size=9,colour = "black"),
        legend.position = "null")+
  geom_signif(y_position = c(11, 11), xmin=c(0.8), xmax=c(1.2),tip_length=0.001,annotations = c("p = 0.001 ***"))+
  ylim(0,12)

##totoal C
TC <- data.input[data.input$ev == "TC",]
bar.TC <- ggplot(data=TC,aes(x=ev, y=mean, fill=water))+
  geom_col(width =0.5, position = "dodge", size=0.05)+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),colour="black",position = position_dodge(0.5), width=0.2,size=0.10)+
  #theme(axis.text.x = element_text(angle))
  labs(x=NULL,y="Soil total C (%)", size=10)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=9, colour = "black"),
        axis.text.y = element_text(size=9,colour = "black"),
        legend.position = "null")+
  geom_signif(y_position = c(1.6, 1.6), xmin=c(0.8), xmax=c(1.2),tip_length=0.001,annotations = c("p = 0.015 **"))+
  ylim(0,1.7)

#biomass
plant.AGB <- data.input[data.input$ev %in% c("AGB","Litter"),]
bar.AGB <- ggplot(data=plant.AGB,aes(x=ev, y=mean, fill=water))+
  geom_col(width =0.5, position = "dodge", size=0.05)+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),colour="black",position = position_dodge(0.5), width=0.2,size=0.10)+
  #theme(axis.text.x = element_text(angle))
  labs(x=NULL,y="Biomass (g m-2)", size=10)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=9, colour = "black"),
        axis.text.y = element_text(size=9,colour = "black"),
        legend.position = "null")+
  geom_signif(y_position = c(7.5, 7.5), xmin=c(0.8), xmax=c(1.2),tip_length=0.001,annotations = c("p = 0.370"))+
  geom_signif(y_position = c(4, 4), xmin=c(1.8), xmax=c(2.2),tip_length=0.001,annotations = c("p = 0.832"))+
  ylim(0,8)

plant.BGB <- data.input[data.input$ev %in% c("SRoot","FRoot"),]
bar.BGB <- ggplot(data=plant.BGB,aes(x=ev, y=mean, fill=water))+
  geom_col(width =0.5, position = "dodge", size=0.05)+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),colour="black",position = position_dodge(0.5), width=0.2,size=0.10)+
  #theme(axis.text.x = element_text(angle))
  labs(x=NULL,y="Biomass (g m-2)", size=10)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=9, colour = "black"),
        axis.text.y = element_text(size=9,colour = "black"),
        legend.position = c(0.8,0.8))+
  geom_signif(y_position = c(0.18, 0.18), xmin=c(0.8), xmax=c(1.2),tip_length=0.001,annotations = c("p = 0.264"))+
  geom_signif(y_position = c(0.11, 0.11), xmin=c(1.8), xmax=c(2.2),tip_length=0.001,annotations = c("p = 0.309"))+
  ylim(0,0.2)
ggarrange(bar.sm,bar.AGB, bar.TC, bar.BGB, nrow = 2, ncol=2, widths = c(1,1.5),align = "hv")

ggarrange(bar.sm, bar.TC, bar.AGB, bar.BGB, nrow = 1, ncol=4, widths = c(1,1,1.5,1.5),align = "h", labels = "auto") %>%
  ggsave(filename = "results/ev.pdf", width = 200, height = 60, units = "mm")

##bar plot of significantly changed genes
library(dplyr)
data.bar <- read.table(file = 'results\\GeoChip_gc1_lnfill0mr_C.stess.virus.gene&split-plot ANOVA_p0.1_for barplot.csv',header=T,sep=",")
head(data.bar)

bar.input <- data.bar %>% 
  subset(select = c(1:3,11:12,17)) %>%
  mutate(pd = cut(P_adjusted, breaks = c(-Inf, 0.001, 0.05, 0.10),
                  labels = c("***", "**", "*")),
         ped = cut(Percentage, breaks = c(-Inf, 0, Inf),
                   labels = c("<0", ">0")),
         sub.category2 = factor(sub.category2, levels = unique(sub.category2)),
         gene.name = factor(gene.name, levels = unique(gene.name)),
         labeld = case_when(Percentage>0 ~ Percentage+se+0.7,
                            Percentage<0 ~ Percentage-se-1.5))
head(bar.input)

range <- bar.input %>% 
  group_by(sub.category2) %>%
  summarise(sub.category1 = unique(sub.category1), start=unique(gene.name)[1], end=last(unique(gene.name))) %>%
  mutate(interval = rep(c("1","0"), length.out = length(start)))

range

input.Cdeg <- bar.input %>%
  subset(sub.category1 == 'Carbon degradation')
range.Cdeg <- range %>%
  subset(sub.category1 == 'Carbon degradation')

##set font family
windowsFonts()
windowsFonts(Times=windowsFont("TT Times New Roman"))
#windowsFonts(Arial=windowsFont("TT Arial"))
#install.packages("extrafont")
#library(extrafont)
#extrafont::font_import(pattern = "arial", prompt = F)
#extrafont::font_import(pattern = "times", prompt = F)
#extrafont::fonttable()
#extrafont::loadfonts(device = "win")

plot.Cdeg <- ggplot(data=input.Cdeg,aes(x=gene.name,y=Percentage, fill=sub.category2))+
  geom_bar(stat = "identity",position="identity",color="black", size=0.05, width = 0.7, alpha = .4)+
  geom_errorbar(aes(ymin=Percentage-se,ymax=Percentage+se),colour="black",width=0.25,size=0.10)+
  scale_fill_brewer(palette="Set2")+
  #theme(axis.text.x = element_text(angle))
  labs(x=NULL,y="Percent change (%)")+
  scale_y_continuous(breaks=seq(-5,25,5), limits=c(-7,27))+
  geom_hline(yintercept=0, size=0.1)+
  geom_text(aes(y=labeld),label=input.Cdeg$pd, size=3, family="Times")+
  theme_bw()+
  theme(legend.position ="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1,vjust=1,size=9, colour = "black"),
        axis.text.y = element_text(size=9,colour = "black"))
plot.Cdeg

##add label of subcategory2
plot.Cdeg.anno1 <- plot.Cdeg +
  geom_text(data = range.Cdeg, aes(y=20, x=end, colour=sub.category2), label=range.Cdeg$sub.category2,
            angle=90, fontface = "bold", size=3)+
  scale_colour_brewer(palette="Set2")
plot.Cdeg.anno1
ggsave(plot.Cdeg.anno1, filename = "results/Cdeg.bar.p0.1.text.pdf", width = 100, height = 80, units = "mm")

library(ggtext)
plot.Cdeg.anno2 <- plot.Cdeg +
  geom_richtext(data = range.Cdeg, aes(y=20, x=end, label=sub.category2), 
                angle=90, fontface = "bold", alpha=.4, size=2)
plot.Cdeg.anno2

ggsave(plot.Cdeg.anno2, filename = "results/Cdeg.bar.p0.1.label.pdf", width = 100, height = 85, units = "mm")

##virus bar plot

input.virus <- bar.input %>%
  subset(sub.category1 == 'Virus')

plot.virus<-ggplot(data=input.virus,aes(x=gene.name,y=Percentage))+
  geom_bar(stat = "identity",position="identity",fill="grey",alpha=0.8,color="black",size=0.05, width = 0.7)+
  geom_errorbar(aes(ymin=Percentage-se,ymax=Percentage+se),colour="black",width=0.25,size=0.10)+
  geom_text(aes(y=labeld),label=input.virus$pd, size=3, family="Times")+
  geom_hline(yintercept=0, size=0.1)+
  labs(x=NULL,y="Percent change (%)")+
  scale_y_continuous(breaks=seq(-5,10,5), limits=c(-7,13))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=9, colour = "black", hjust = 1, vjust = 1, angle = 30),
        axis.text.y = element_text(size=9,colour = "black"))+
  facet_grid(. ~ sub.category2, scale="free_x",space="free")
plot.virus
ggsave(plot.virus, filename = "results/virus bar_lnmr_1.pdf", width=70,height=75, units = "mm")


###bar plot of CH4
##C degradation bar plot
input.methane <- bar.input %>%
  subset(sub.category1 == 'Methane') %>%
  arrange(desc(sub.category2)) %>%
  mutate(sub.category2=factor(sub.category2, levels=unique(sub.category2))) %>%
  mutate(gene.name = factor(gene.name, levels=unique(gene.name)))

  
CH4bar<-ggplot(data=input.methane,aes(x=gene.name,y=Percentage,fill=sub.category2))+
  geom_bar(stat = "identity",position="identity", size=0.05, width = 0.7)+
  geom_errorbar(aes(ymin=Percentage-se,ymax=Percentage+se),colour="black",width=0.25,size=0.10)+
  scale_fill_brewer(palette="Pastel2")+
  geom_text(aes(y=labeld),label=input.methane$pd, size=3, family="Times")+
  labs(x=NULL,y="Percent change (%)", fill="Subcategory")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,vjust=1,size=9, colour = "black"),
        axis.text.y = element_text(size=9,colour = "black"))
CH4bar

ggsave(CH4bar, filename = "results/CH4 bar_lnmr.pdf", width=4,height=3)

library(ggpubr)
ggarrange(plot.Cdeg.anno1, plot.virus,nrow = 1, ncol=2, widths = c(1.8,1), labels = c("e","f")) %>%
  ggsave(filename = "results/Cdeg&virus.lnmr2.pdf", width = 200, height = 80, units = "mm")


##bar plot of stress
input.stress <- bar.input %>%
  subset(sub.category1 == 'Stress')
range.stress <- range %>%
  subset(sub.category1 == 'Stress')

windowsFonts()
windowsFonts(Times=windowsFont("TT Times New Roman"))

plot.stress <- ggplot(data=input.stress,aes(x=gene.name,y=Percentage, fill=sub.category2))+
  geom_bar(stat = "identity",position="identity",color="black", size=0.05, width = 0.7, alpha = .4)+
  geom_errorbar(aes(ymin=Percentage-se,ymax=Percentage+se),colour="black",width=0.25,size=0.10)+
  scale_fill_brewer(palette="Set2")+
  #theme(axis.text.x = element_text(angle))
  labs(x=NULL,y="Percent change (%)", fill="Subcategory")+
  #scale_y_continuous(limits = c(0,22))+
  geom_hline(yintercept=0, size=0.1)+
  geom_text(aes(y=labeld),label=input.stress$pd, size=3, family="Times")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1,vjust=1,size=9, colour = "black"),
        axis.text.y = element_text(size=9,colour = "black"))

ggsave(plot.stress, filename = "results/stress.bar_lnmr_p0.1_legend.pdf", width = 180, height =70, units = "mm")

  
 plot.stress.text <- ggplot(data=input.stress,aes(x=gene.name,y=Percentage, fill=sub.category2))+
    geom_bar(stat = "identity",position="identity",color="black", size=0.05, width = 0.7, alpha = .4)+
    geom_errorbar(aes(ymin=Percentage-se,ymax=Percentage+se),colour="black",width=0.25,size=0.10)+
    scale_fill_brewer(palette="Set2")+
    #theme(axis.text.x = element_text(angle))
    labs(x=NULL,y="Percent change (%)")+
    scale_y_continuous(limits = c(0,33))+
    geom_hline(yintercept=0, size=0.1)+
    geom_text(aes(y=labeld),label=input.stress$pd, size=3, family="Times")+
    theme_bw()+
    theme(legend.position ="none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1,vjust=1,size=9, colour = "black"),
          axis.text.y = element_text(size=9,colour = "black"))+
    geom_text(data = range.stress, aes(y=27, x=end, colour=sub.category2), label=range.stress$sub.category2,
              angle=90, fontface = "bold", size=3)+
    scale_colour_brewer(palette="Set2")
ggsave(plot.stress.text, filename = "results/stress.bar.p0.1.text.pdf", width = 120, height =100, units = "mm")




##bar plot of bacteriophage lysis gene
data.lysis<-read.table('results\\GeoChip_gc1_lnfill0mr_phage lysis.gene&split-plot ANOVA_for ggplot.csv', header = T, sep = ",")
head(data.lysis)
input.lysis<- data.lysis %>%
  mutate(pd = cut(p.value, breaks = c(-Inf, 0.001, 0.05, 0.10),
                  labels = c("***", "**", "*")),
         ped = cut(Percentage, breaks = c(-Inf, 0, Inf),
                   labels = c("<0", ">0")),
          gene.name = factor(gene.name, levels = unique(gene.name)),
         labeld = case_when(Percentage>0 ~ Percentage+se+1.5,
                            Percentage<0 ~ Percentage-se-1.5))

plot.lysis<-ggplot(data=input.lysis,aes(x=gene.name,y=Percentage))+
  geom_bar(stat = "identity",position="identity",fill="grey",alpha=0.8,color="black",size=0.05, width = 0.7)+
  geom_errorbar(aes(ymin=Percentage-se,ymax=Percentage+se),colour="black",width=0.25,size=0.10)+
  geom_text(aes(y=labeld),label=input.lysis$pd, size=3, family="Times")+
  geom_hline(yintercept=0, size=0.1)+
  labs(x=NULL,y="Percent change (%)")+
  scale_y_continuous(breaks=seq(-15,30,10), limits=c(-18,30))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=8.5, colour = "black", hjust = 1, vjust = 1, angle = 45),
        axis.text.y = element_text(size=9,colour = "black"))
plot.lysis
ggsave(plot.lysis, filename = "results/bacteriophage lysis gene_lnmr_barplot.pdf", width=70,height=80, units = "mm")
