library(data.table)
library(DESeq2)
library(ggplot2)
library(ggrepel)
setwd('/Users/malleyce/Downloads/DRGs/')
load('/Users/malleyce/Downloads/DRGs/DRGs.Nociceptors.DDS.RData')

#how to subset
#dds.subset <- dds[,dds@colData@listData$condition %in% c('H9noc_D28','GTEx_Hypothalamus','GTEx_Hippocampus','GTEx_Cortex',"GTEx_Spinal_cord_cervical_c-1")]
dds.subset <- dds

matrixFile <- as.data.frame(counts(dds.subset, normalized=FALSE))
mat <- as.matrix(matrixFile)

lowCountLimit <- 25
minSamplesLimit <- 3
factorName <- "condition"

normCounts <- data.frame(counts(dds.subset, normalized=TRUE))

normCounts[,"goodCount"]<-rowSums(normCounts > lowCountLimit)

redMat <- subset(normCounts, goodCount >= minSamplesLimit)
redMat <- redMat[,1:(ncol(redMat)-1)]
conds <- factor(dds.subset@colData$condition)

pca <- prcomp(t(redMat), center=TRUE, scale=TRUE)
conditions <- data.frame("Condition"= conds)

pca.df <- as.data.frame(pca$x)

sampleTable <- as.data.table(readxl::read_xlsx('/Users/malleyce/Downloads/DRGs/DRGs_merged_sampletable.xlsx'))

ggplot(pca.df, aes(x=pca$x[,"PC1"], y=pca$x[,"PC2"])) + geom_point(aes(color=conditions$Condition), size=5, alpha=0.5) + theme_bw() +
  geom_text_repel(aes(label=sampleTable$sample_id),color="black") +
  labs(x="PC1", y="PC2", title="PCA of DRGs plus Tao's nociceptors") +
  guides(color=guide_legend(title="Condition"))+theme(panel.grid=element_line(size=1), axis.text = element_text(size=12), axis.title=element_text(face='plain'), title=element_text(face = 'bold'), legend.text = element_text(size=12), legend.title = element_text(face='plain', size=12))

# merge sample metadata-----
# 
# sratable <- as.data.table(read_xlsx('SraRunaTablesubset.xlsx'))
# sratable
# 
# sample_multi <- fread('dbgap_pheno/phs001158.v2.pht005593.v2.p1.Human_Dorsal_Root_Ganglion_Sample.MULTI.txt', skip=10L)
# subject_multi <- fread('dbgap_pheno/phs001158.v2.pht005592.v2.p1.Human_Dorsal_Root_Ganglion_Subject.MULTI.txt', skip=10L)
# 
# dbgap_info <- merge(sample_multi, subject_multi, by=c('dbGaP_Subject_ID', 'SUBJECT_ID'))
# dbgap_info
# names(dbgap_info)[4] <- 'BioSample'
# 
# sratable_merged <- merge(dbgap_info, sratable, by=c('BioSample'))
# sratable_merged
# 
# exclude <-c('SRR8533972','SRR8533976','SRR8533978')
# 
# cat(row.names(as.data.frame(dds@colData)), sep='\n')

# add other metadata-----


coldata <- as.data.frame(dds@colData)
coldata$sample <- row.names(coldata)
coldata <- as.data.table(coldata)
coldata$order <- 1:nrow(coldata)

sampleTable <- as.data.table(readxl::read_xlsx('/Users/malleyce/Downloads/DRGs/DRGs_merged_sampletable.xlsx'))
sampletable <- sampleTable

coldata <- merge(coldata, sampletable, by.x=c('sample', 'condition'), by.y=c('sampleFile', 'condition'))
coldata <- coldata[order(order)]
coldata

plottable <- pca.df
plottable$file <- coldata$fileName
plottable <- as.data.table(plottable)
plottable <- plottable[,c('file','PC1','PC2')]
plottable$condition <- coldata$condition
plottable$source <- coldata$sample_source
plottable$replicate <- coldata$sample_id
plottable$liblayout <- coldata$library_layout
plottable$patient_sex <- coldata$patient_sex
plottable$condition_relevel <- factor(plottable$condition, levels=c('H9noc_D0', 'H9noc_D4', 'H9noc_D8', 'H9noc_D12', 'H9noc_D21', 'H9noc_D28', 'H9noc_D56', 'DRG_chronic_pain', 'DRG_control'))

attach(plottable)

pca.plot <- ggplot(plottable, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=condition_relevel,
                 text = paste(
                   "Sample: ", condition_relevel, "\n",
                   "Replicate: ", replicate, "\n",
                   "Source: ", source, "\n",
                   "Patient sex: ", patient_sex, "\n",
                   "Library layout: ", liblayout, "\n",
                   sep = ""
                 )), size=5, alpha=0.5) + 
  theme_bw() +
  labs(x="PC1 (31%)", y="PC2 (21%)", 
       title="PCA of nociceptor samples plus dbGaP DRGs") +
  guides(color=guide_legend(title="Sample"))+
  theme(panel.grid=element_line(size=1),
        axis.text = element_text(size=12), axis.title=element_text(face='plain'),
        title=element_text(face = 'bold'), legend.text = element_text(size=12),
        legend.title = element_text(face='plain', size=12))

library("plotly")
library("tidyverse")
library("htmlwidgets")

ggplotly(pca.plot, tooltip = "text")
#