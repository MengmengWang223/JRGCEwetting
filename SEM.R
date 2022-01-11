library(lavaan)

data <- read.csv("ENV and comm data for SEM.csv", row.names = 1, header = T)
dim(data)
root <- data$FRoot + data$SRoot
data <- cbind(data, root)
## scale to the same variance
dat <- data.frame(scale(data, center = F))
pairs(dat)
colnames(dat)
##only at AN
data1 <- rbind(data[1:16,], data[33:48,])
dat1 <- data.frame(scale(data1, center = F))
##data at EN
data2 <- rbind(data[17:32,], data[49:64,])
dat2 <- data.frame(scale(data2, center = F))


# mode3, add water treatment
model3 <- '
# regressions
Moisture ~ water
pH ~ water
FRoot ~  pH
Litter ~ Moisture
TC ~ Moisture + FRoot+ sigchg.recalC+PCo1.OTU
sigchg.recalC ~ phage + water
PCo1.OTU ~  pH + Litter
phage ~ water + Litter + PCo1.GeoChip.ln0mr
PCo1.GeoChip.ln0mr ~ diversity.GeoChip + Moisture
diversity.GeoChip ~  FRoot + pH + AGB
AGB ~  AG
AG ~  water + Moisture
'

carbonsem3 <- sem(model3, missing = "ml", data = dat1)

summary(carbonsem3, rsquare = T, standardized = T, fit = T)

residuals(carbonsem3, type = "cor")
modificationIndices(carbonsem3, standardized = F)
