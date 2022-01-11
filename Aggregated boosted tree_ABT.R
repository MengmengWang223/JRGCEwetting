##Aggregated boosted tree, ABT
#gbm package
library(gbm)
library(picante)
library(ape)
library(vegan)
library(ade4)
library(caret)
library(ggplot)

##----------------------------------------------------------------------------
##ABT for soil C
##input prepared data
##ABT for soil carbon
data <- read.csv("ENV and comm data for SEM.csv", row.names = 1, header = T)
dim(data)
head(data)
data.1 <- data[,colnames(data) %in% c("PCo1.OTU", "Litter", "AGB","PCo1.GeoChip.ln0mr","diversity.OTU", "sigchg.recalC", "diversity.GeoChip", "phage")]
colnames(data.1)

##AN
data.AN <- data.1[c(1:16,33:48),]
dim(data.AN)
head(data.AN)

## scale to the same variance
dat <- data.frame(scale(data.AN, center = F))
dim(dat)
pairs(dat)
colnames(dat)

##For soil C
gbm.C <- gbm(TC ~., data=dat, distribution="gaussian", n.trees=10000, shrinkage=0.001,
             bag.fraction = 0.85, cv.folds=5, interaction.depth = 1, verbose=FALSE)

# get MSE and compute RMSE
sqrt(min(gbm.C$cv.error))
#Check performance using the cv
best.iter.C <- gbm.perf(gbm.C, method="cv")
print(best.iter.C)

# Plot relative influence of each variable
par(mar = c(5, 9, 1, 1))
summary(gbm.C, n.trees = best.iter.C,las=1) # using estimated best number of trees

##create hyperparameter grid
hyper_grid <- expand.grid(
  shrinkage = c(.001, .01),
  interaction.depth = c(1, 3, 5),
  n.minobsinnode = c(3, 5, 10),
  bag.fraction = c(0.85, 0.9),
  optimal_trees = 0,               # a place to dump results
  min_RMSE = 0                     # a place to dump results
)

# total number of combinations
nrow(hyper_grid)
#loop grid research
for(i in 1:nrow(hyper_grid)) {
  # train model
  gbm.tune <- gbm(
    TC ~ ., data=dat,
    distribution = "gaussian",
    n.trees = 30000,
    interaction.depth = hyper_grid$interaction.depth[i],
    shrinkage = hyper_grid$shrinkage[i],
    n.minobsinnode = hyper_grid$n.minobsinnode[i],
    bag.fraction = hyper_grid$bag.fraction[i],
    cv.folds = 5,
    n.cores = NULL, # will use all cores by default
    verbose = FALSE
  )
  
  # add min training error and trees to grid
  hyper_grid$optimal_trees[i] <- which.min(gbm.tune$cv.error)
  hyper_grid$min_RMSE[i] <- sqrt(min(gbm.tune$cv.error))
}

hyper_grid %>%
  dplyr::arrange(min_RMSE) %>%
  head(20)

write.csv(hyper_grid, file="results/ABT_soilC/hyper_grid.csv")

gbm.C.final <- gbm(TC ~., data=dat, distribution="gaussian", n.trees=50000, shrinkage=0.01,
             bag.fraction = 0.85, cv.folds=5, interaction.depth = 1, n.minobsinnode =3, verbose=FALSE)

sqrt(min(gbm.C.final$cv.error))
best.iter.C.final <- gbm.perf(gbm.C.final, method="cv")
print(best.iter.C.final)

# Plot relative influence of each variable
par(mar = c(5, 10, 1, 1))
summary(
  gbm.C.final,
  method = relative.influence, # also can use permutation.test.gbm
  las = 1
)

#model evaluation
fitControl=trainControl(method="cv", number=10, returnResamp="all")
model3=train(TC~., data=dat, distribution="gaussian", method='gbm', bag.fraction=0.85, trControl=fitControl, 
             verbose=F, tuneGrid=data.frame(n.trees=best.iter.C.final, shrinkage=0.01, interaction.depth=1, n.minobsinnode=3))
model3

##plot
ABT.TC.in <- read.delim("clipboard")

ABT.TC.in <- ABT.TC.in %>%
  arrange(desc(rel.inf)) %>%
  mutate(var=factor(var, levels=unique(var)))
dim(ABT.TC.in)
head(ABT.TC.in)

windowsFonts()
windowsFonts(Times=windowsFont("TT Times New Roman"))
ev.name <- c("Fine root", "Moisture", "Recal-C degra. gene", "pH", "Taxonomic composition", "Litter",
             "Bacteriophage gene", "AGB", "Functional composition", "Shallow root")
plot.ABT.TC.in <- ggplot(data=ABT.TC.in, aes(x=var, y=rel.inf))+
  geom_bar(stat = "identity",fill="#beaed4", alpha=.8, width = .7)+
  geom_text(aes(y=rel.inf+1.5),label=round(ABT.TC.in$rel.inf,digits = 1), size=3, family="Times")+
  #coord_flip()+
  labs(x=NULL,y="Relative influence (%)")+
  scale_y_continuous(limits = c(0,24), expand = c(0, 0))+
  scale_x_discrete(labels=ev.name)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,vjust=1,size=9, colour = "black"),
        axis.text.y = element_text(size=9,colour = "black"))

ggsave(plot.ABT.TC.in, filename = "results/ABT_soilC/ABT_C_ggplot_vertical.pdf", width = 70, height = 85, units = "mm")