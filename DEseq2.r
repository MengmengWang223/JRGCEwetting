## RNA-seq analysis with DESeq2
## Stephen Turner, @genetics_blog

# RNA-seq data from GSE52202
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse52202. All patients with
# ALS, 4 with C9 expansion ("exp"), 4 controls without expansion ("ctl")

# Import & pre-process ----------------------------------------------------

# Import data from featureCounts
## Previously ran at command line something like this:
## featureCounts -a genes.gtf -o counts.txt -T 12 -t exon -g gene_id GSM*.sam
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

setwd("D:\\lab_Tsinghua\\aProjects\\aStanford\\water_32 vs 32")

OTU <- read.table("OTU_nosinglet_morethan0.001%.csv", header = TRUE, sep = ",", row.names = 1)
dim(OTU)
OTU.anno <- OTU[,1:10]
OTU.data <- OTU[,-(1:10)]
dim(OTU.data)

countdata <- OTU.data
head(countdata)
dim(countdata)


# Assign condition (first four are controls, second four contain the expansion)
(condition <- factor(c(rep("A", 32), rep("B", 32))))

# Analysis with DESeq2 ----------------------------------------------------
library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata1 <- data.frame(row.names=colnames(countdata), condition))

ddsdata <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata1, design=~condition)
ddsdata

# Run the DESeq pipeline
dds <- DESeq(ddsdata)


# Get differential expression results
res <- results(dds)
table(res$padj<0.05)

## Order by adjusted p-value
res <- res[order(res$padj), ]

## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=FALSE)), by="row.names", sort=FALSE)
head(resdata)
resdata2 <- resdata[,-1]
rownames(resdata2) <- resdata[,1]
head(resdata2)
resdata.anno <- merge(resdata2, OTU.anno, by="row.names", sort=FALSE)
head(resdata.anno)

## Write results
write.csv(resdata.anno, file="results\\OTU&DESeq2&split-plot ANOVA\\OTU_morethan0.001%_DESeq2 results.csv")
