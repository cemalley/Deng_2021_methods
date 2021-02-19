library(DESeq2)
library(data.table)
library(RUVSeq) # http://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf
library(EDASeq)
library(readxl)

#setwd("/data/NCATS_ifx/data/dbGaP/DRGs/SRA/DESeq2")
setwd("~/Downloads/DRGs/")

sampleTable <- as.data.table(readxl::read_xlsx('DRGs_merged_sampletable.xlsx'))

directory <- "~/Downloads/DRGs/Countfiles_gene_symbol/intersection/"
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = ~ condition)

keep <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)
setwd("~/Downloads/DRGs/")
save(dds, file='DRGs.Nociceptors.DDS.RData')

load('DRGs.Nociceptors.DDS.RData')

x <- as.factor(sampleTable$condition)
set <- newSeqExpressionSet(as.matrix(counts(dds, normalized=F)),
                           phenoData = data.frame(x, row.names= c(sampleTable$sample_id)))

spikes <- unique(c('ANAPC5', 'ANAPC15', 'ARID3B', 'ARL10', 'ATXN2', 'C16orf62', 'C3orf49', 'CCAR1', 'CCDC125', 'CCDC90B', 'CHFR', 'DHRSX', 'FRMD8', 'GGA1', 'HERC4', 'MKNK1', 'NASP', 'NME4', 'OTUB1', 'PMF1', 'POLR2B', 'POLR3A', 'POMK', 'PSMA3-AS1', 'PTPN14', 'RAPGEF6', 'REL', 'RRP1', 'RUNDC1', 'SAMD4B', 'SLC4A1AP', 'SLMAP', 'SMARCAL1', 'SNAP29', 'SNRNP200', 'SUPT4H1', 'TBC1D22A', 'THUMPD3-AS1', 'TSPOAP1-AS1', 'TUBGCP2', 'WDTC1', 'ZNF544','C1orf43','CHMP2A','EMC7','GPI','PSMB2','PSMB4','RAB7A','REEP5','SNRPD3','VCP','VPS29')) #https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-019-0538-z#Tab3

spikes <- spikes[spikes %in% row.names(dds)]

set2 <- RUVg(set, spikes, k=1)

condition <- x

names(pData(set2)) <- c('condition', 'W_1')

dds <- DESeqDataSetFromMatrix(countData = counts(set2),
                              colData = pData(set2),
                              design = ~ W_1 + condition)
dds <- DESeq(dds)

#setwd("/data/NCATS_ifx/data/dbGaP/DRGs/SRA/DESeq2")
save(dds, file='DRGs.Nociceptors.DDS.RData')


quit(save='no')

# sbatch --mem=5g --time=10-0:00:00 make_DDS.sh #  did not take more than 3G ram in this set.
# in /data/NCATS_ifx/data/dbGaP/DRGs/SRA/DESeq2

### redoing to add day 56 and rerun ruvseq-----
load('~/Downloads/Noci_GTEx_DRGs/Noci_GTEx_DRGs_corrected_DDS.RData')
dds.noci <- dds[ , dds$condition %in% c('DRG_chronic_pain', 'DRG_control','H9noc_D0', 'H9noc_D4', 'H9noc_D8', 'H9noc_D12', 'H9noc_D21', 'H9noc_D28', 'H9noc_D56')]
cts.raw <- as.data.frame(counts(dds.noci, normalized=F))
cts.raw
dds.noci@colData

conditions <- as.data.frame(dds.noci@colData)
conditions$relevel <- factor(conditions$condition, levels=c('DRG_chronic_pain', 'DRG_control','H9noc_D0', 'H9noc_D4', 'H9noc_D8', 'H9noc_D12', 'H9noc_D21', 'H9noc_D28', 'H9noc_D56'))

x <- as.factor(conditions$relevel)
set <- newSeqExpressionSet(as.matrix(counts(dds.noci, normalized=F)),
                           phenoData = data.frame(x, row.names= row.names(dds.noci@colData)))

spikes <- unique(c('ANAPC5', 'ANAPC15', 'ARID3B', 'ARL10', 'ATXN2', 'C16orf62', 'C3orf49', 'CCAR1', 'CCDC125', 'CCDC90B', 'CHFR', 'DHRSX', 'FRMD8', 'GGA1', 'HERC4', 'MKNK1', 'NASP', 'NME4', 'OTUB1', 'PMF1', 'POLR2B', 'POLR3A', 'POMK', 'PSMA3-AS1', 'PTPN14', 'RAPGEF6', 'REL', 'RRP1', 'RUNDC1', 'SAMD4B', 'SLC4A1AP', 'SLMAP', 'SMARCAL1', 'SNAP29', 'SNRNP200', 'SUPT4H1', 'TBC1D22A', 'THUMPD3-AS1', 'TSPOAP1-AS1', 'TUBGCP2', 'WDTC1', 'ZNF544','C1orf43','CHMP2A','EMC7','GPI','PSMB2','PSMB4','RAB7A','REEP5','SNRPD3','VCP','VPS29')) #https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-019-0538-z#Tab3

spikes <- spikes[spikes %in% row.names(dds)]

set2 <- RUVg(set, spikes, k=1)

condition <- x

names(pData(set2)) <- c('condition', 'W_1')

counts.set2 <- as.data.frame(counts(set2))
counts.set2 <- subset(counts.set2, select=c('Sample_1_H9noc_D0_1_intersection_gene_symbol_counts.txt', 'Sample_2_H9noc_D0_2_intersection_gene_symbol_counts.txt', 'Sample_3_H9noc_D0_3_intersection_gene_symbol_counts.txt', 'Sample_4_H9noc_D4_1_intersection_gene_symbol_counts.txt', 'Sample_5_H9noc_D4_2_intersection_gene_symbol_counts.txt', 'Sample_6_H9noc_D4_3_intersection_gene_symbol_counts.txt', 'Sample_7_H9noc_D8_1_intersection_gene_symbol_counts.txt', 'Sample_8_H9noc_D8_2_intersection_gene_symbol_counts.txt', 'Sample_9_H9noc_D8_3_intersection_gene_symbol_counts.txt', 'Sample_10_H9noc_D12_1_intersection_gene_symbol_counts.txt', 'Sample_11_H9noc_D12_2_intersection_gene_symbol_counts.txt', 'Sample_12_H9noc_D12_3_intersection_gene_symbol_counts.txt', 'Sample_13_H9noc_D21_1_intersection_gene_symbol_counts.txt', 'Sample_14_H9noc_D21_2_intersection_gene_symbol_counts.txt', 'Sample_15_H9noc_D21_3_intersection_gene_symbol_counts.txt', 'H9_Noc0430_D28_1_intersection_gene_symbol_counts.txt', 'H9_Noc0430_D28_2_intersection_gene_symbol_counts.txt', 'H9_Noc0430_D28_3_intersection_gene_symbol_counts.txt', 'H9_Noc_D56_1_intersection_gene_symbol_counts.txt', 'H9_Noc_D56_2_intersection_gene_symbol_counts.txt', 'H9_Noc_D56_3_intersection_gene_symbol_counts.txt', 'SRR8533960_intersection_gene_symbol_counts.txt', 'SRR8533961_intersection_gene_symbol_counts.txt', 'SRR8533962_intersection_gene_symbol_counts.txt', 'SRR8533963_intersection_gene_symbol_counts.txt', 'SRR8533964_intersection_gene_symbol_counts.txt', 'SRR8533965_intersection_gene_symbol_counts.txt', 'SRR8533966_intersection_gene_symbol_counts.txt', 'SRR8533967_intersection_gene_symbol_counts.txt', 'SRR8533968_intersection_gene_symbol_counts.txt', 'SRR8533969_intersection_gene_symbol_counts.txt', 'SRR8533970_intersection_gene_symbol_counts.txt', 'SRR8533971_intersection_gene_symbol_counts.txt', 'SRR8533973_intersection_gene_symbol_counts.txt', 'SRR8533974_intersection_gene_symbol_counts.txt', 'SRR8533975_intersection_gene_symbol_counts.txt', 'SRR8533977_intersection_gene_symbol_counts.txt', 'SRR8533979_intersection_gene_symbol_counts.txt', 'SRR8533980_intersection_gene_symbol_counts.txt', 'SRR8533981_intersection_gene_symbol_counts.txt', 'SRR8533982_intersection_gene_symbol_counts.txt', 'SRR8533983_intersection_gene_symbol_counts.txt', 'SRR8533984_intersection_gene_symbol_counts.txt', 'SRR8533985_intersection_gene_symbol_counts.txt', 'SRR8533986_intersection_gene_symbol_counts.txt'))
names(counts.set2)

dds <- DESeqDataSetFromMatrix(countData = counts.set2,
                              colData = colData,
                              design = ~ W_1 + condition)
dds <- DESeq(dds)
#colData <- as.data.frame(dds@colData)
#write.table(colData, '/Users/malleyce/Downloads/DRGs/coldata.temp.csv')
#colData <- read.csv('/Users/malleyce/Downloads/DRGs/coldata.temp.csv', row.names=1)
#dds@colData <- colData
View(as.data.frame(dds@colData))

save(dds, file='/Users/malleyce/Downloads/DRGs/DRGs.Nociceptors.DDS.RData')

# DE tests-----
library(doParallel)

source('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/DESeq2-pipeline-function-prototype.R')

#tests_table <- data.table(condition1=c('H9noc_D56','H9noc_D56'), condition2=c('DRG_chronic_pain','DRG_control'))
tests_table <- data.table(condition1=c('H9noc_D28','H9noc_D28'), condition2=c('DRG_chronic_pain','DRG_control'))

testdir <- '/Users/malleyce/Downloads/DRGs/DE'
volcanodir <- '/Users/malleyce/Downloads/DRGs/DE/Volcano'

keep <- rowSums( counts(dds, normalized=TRUE) >= 20 ) >= 3
dds <- dds[keep,]

cl <- makeCluster(2)
registerDoParallel(cl)
foreach(i=1:nrow(tests_table), .packages='data.table') %dopar% {
  condition1 <- tests_table$condition1[i]
  condition2 <- tests_table$condition2[i]
  
  testdir <- testdir
  volcanodir <- volcanodir
  
  DESeq2_pipeline(dds, condition1, condition2, testdir, volcanodir)
}

stopCluster(cl)

# volcano plot, noci D56 versus DRGs pain-----
library(data.table)
library(ggthemes)
library(ggrepel)
library(cowplot)

res <- fread('/Users/malleyce/Downloads/DRGs/DE/DE.H9noc_D56.vs.DRG_chronic_pain.csv')
restoplot <- na.omit(res)
condition1 <- 'H9noc_D56'
condition2 <- 'DRG_chronic_pain'

# set p to 1*10^-300 where padj=0 for plotting gene labels more easily.
restoplot[,padj_reduce:=padj]
restoplot[padj_reduce ==0,padj_reduce := (1*10^(-300))]

right <- subset(restoplot, ( ((log2FoldChange > 1) & (-log10(padj_reduce) > 10)))  ) #0
left <- subset(restoplot, ( ((log2FoldChange < -1) & (-log10(padj_reduce) > 10)))) #0
down.nrow <- nrow(left)
up.nrow <- nrow(right)

if(nrow(right)>15){
  right <- right[order(-log2FoldChange)]
  right <- right[1:15,]
}
#right <- right[GeneId %in% c('C1QC', 'TINAGL1', 'C1QB', 'C1QA', 'STAB1', 'IGHG1', 'GJA4', 'TYROBP', 'STEAP4', 'RGS1', 'TSPAN8', 'MEG3', 'VSIG4', 'CYP4B1', 'FCGR3A','HTR3A', 'C3', 'NTRK1', 'RARRES2', 'GPX2', 'TRARG1', 'TMEM176A', 'SFRP5', 'NCMAP', 'IL31RA', 'SERPING1', 'SLC15A3', 'RASL12', 'NEFH', 'PPM1J', 'STEAP3', 'SPP1', 'AVIL', 'CPT1A', 'OSMR', 'AIF1L', 'C2CD2'),]

if(nrow(left)>15){
  left <- left[order(log2FoldChange)]
  left <- left[1:15,]
}

restoplot[,threshold:=ifelse( ((-log10(padj_reduce) > 10) & log2FoldChange < -1),
                             '#007F00', 
       ifelse( ((-log10(padj_reduce) > 10) & log2FoldChange > 1),
               'red', 'gray' ))]

plot <- ggplot(data=restoplot) + geom_point(aes(x=log2FoldChange, y=-log10(padj_reduce)),
                                            color=restoplot$threshold)+
  geom_text_repel(data=right, aes(x=log2FoldChange, y=-log10(padj_reduce), label=GeneId),
                  segment.size = 0.5)+
  geom_text_repel(data=left, aes(x=log2FoldChange, y=-log10(padj_reduce), label=GeneId),
                  segment.size = 0.5)+
  geom_vline(aes(xintercept=1),color='brown')+
  geom_vline(aes(xintercept=-1),color='brown')+
  geom_hline(aes(yintercept=10), color='brown')+
  labs(title=paste0('Volcano plot: ', condition1, ' vs ' ,
                    condition2, '. UP: ',up.nrow,', DOWN: ',down.nrow), x='Effect size: log2(fold-change)',
       y='-log10(adjusted p-value)')+
  theme_bw()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12))+
  geom_text(label="Higher in DRG chronic pain", x=15, y=110, color='red', vjust=1)+
  geom_text(label="Higher in H9noc D56", x=-16, y=110, color='#007F00', vjust=1)+
  xlim(-20,20)+ylim(0,120) # may need to expand to -20,20 if the log2fold changes are that big
plot

# volcano plot, noci D56 versus DRGs control-----
library(data.table)
library(ggthemes)
library(ggrepel)
library(cowplot)

res <- fread('/Users/malleyce/Downloads/DRGs/DE/DE.H9noc_D56.vs.DRG_control.csv')
restoplot <- na.omit(res)
condition1 <- 'H9noc_D56'
condition2 <- 'DRG_control'

# set p to 1*10^-300 where padj=0 for plotting gene labels more easily.
restoplot[,padj_reduce:=padj]
restoplot[padj_reduce ==0,padj_reduce := (1*10^(-300))]

right <- subset(restoplot, ( ((log2FoldChange > 1) & (-log10(padj_reduce) > 10)))  ) #0
left <- subset(restoplot, ( ((log2FoldChange < -1) & (-log10(padj_reduce) > 10)))) #0
down.nrow <- nrow(left)
up.nrow <- nrow(right)

if(nrow(right)>15){
  right <- right[order(-log2FoldChange)]
  right <- right[1:15,]
}

if(nrow(left)>15){
  left <- left[order(log2FoldChange)]
  left <- left[1:15,]
}

restoplot[,threshold:=ifelse( ((-log10(padj_reduce) > 10) & log2FoldChange < -1),
                              '#007F00', 
                              ifelse( ((-log10(padj_reduce) > 10) & log2FoldChange > 1),
                                      'red', 'gray' ))]

plot <- ggplot(data=restoplot) + geom_point(aes(x=log2FoldChange, y=-log10(padj_reduce)),
                                            color=restoplot$threshold)+
  geom_text_repel(data=right, aes(x=log2FoldChange, y=-log10(padj_reduce), label=GeneId),
                  segment.size = 0.5)+
  geom_text_repel(data=left, aes(x=log2FoldChange, y=-log10(padj_reduce), label=GeneId),
                  segment.size = 0.5)+
  geom_vline(aes(xintercept=1),color='brown')+
  geom_vline(aes(xintercept=-1),color='brown')+
  geom_hline(aes(yintercept=10), color='brown')+
  labs(title=paste0('Volcano plot: ', condition1, ' vs ' ,
                    condition2, '. UP: ',up.nrow,', DOWN: ',down.nrow), x='Effect size: log2(fold-change)',
       y='-log10(adjusted p-value)')+
  theme_bw()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12))+
  geom_text(label="Higher in DRG control", x=15, y=120, color='red', vjust=1)+
  geom_text(label="Higher in H9noc D56", x=-16, y=120, color='#007F00', vjust=1)+
  xlim(-20,20)+ylim(0,120) # may need to expand to -20,20 if the log2fold changes are that big
plot

# what is the difference of DE genes between chronic and control-----
res.DRGs <- fread('/Users/malleyce/Downloads/DRGs/DE/DE.H9noc_D56.vs.DRG_chronic_pain.csv')
res.control <- fread('/Users/malleyce/Downloads/DRGs/DE/DE.H9noc_D56.vs.DRG_control.csv')

all.genes <- unique(c(res.DRGs$GeneId, res.control$GeneId))
library(Hmisc)

names(res.DRGs) <- c('GeneId','H9noc_D56','DRG_chronic_pain',
                     'L2FC_DRGcp','lfcSE_DRGcp','stat_DRGcp',
                     'pval_DRGcp','padj_DRGcp')

names(res.control) <- c('GeneId','H9noc_D56','DRG_control',
                     'L2FC_DRGct','lfcSE_DRGct','stat_DRGct',
                     'pval_DRGct','padj_DRGct')


res.merged <- merge(res.DRGs, res.control, by=c('GeneId','H9noc_D56'))
res.merged <- res.merged[,c('GeneId','H9noc_D56','DRG_chronic_pain',
                            'DRG_control','L2FC_DRGcp','L2FC_DRGct',
                            'padj_DRGcp','padj_DRGct')]

res.merged[,padj_diff:=padj_DRGcp-padj_DRGct]
res.merged <- res.merged[order(padj_diff)]
cat(res.merged$GeneId[1:100], sep='\n')
fwrite(res.merged, 'Difference_DE_H9nociD56_DRGs.csv')


# volcano plot, noci D28 versus DRGs pain-----
library(data.table)
library(ggthemes)
library(ggrepel)
library(cowplot)

res <- fread('/Users/malleyce/Downloads/DRGs/DE/DE.H9noc_D28.vs.DRG_chronic_pain.csv')
restoplot <- na.omit(res)
condition1 <- 'H9noc_D28'
condition2 <- 'DRG_chronic_pain'

# set p to 1*10^-300 where padj=0 for plotting gene labels more easily.
restoplot[,padj_reduce:=padj]
restoplot[padj_reduce ==0,padj_reduce := (1*10^(-300))]

right <- subset(restoplot, ( ((log2FoldChange > 1) & (-log10(padj_reduce) > 10)))  ) #0
left <- subset(restoplot, ( ((log2FoldChange < -1) & (-log10(padj_reduce) > 10)))) #0
down.nrow <- nrow(left)
up.nrow <- nrow(right)

if(nrow(right)>15){
  right <- right[order(-log2FoldChange)]
  right <- right[1:15,]
}

if(nrow(left)>15){
  left <- left[order(log2FoldChange)]
  left <- left[1:15,]
}

restoplot[,threshold:=ifelse( ((-log10(padj_reduce) > 10) & log2FoldChange < -1),
                              '#007F00', 
                              ifelse( ((-log10(padj_reduce) > 10) & log2FoldChange > 1),
                                      'red', 'gray' ))]

plot <- ggplot(data=restoplot) + geom_point(aes(x=log2FoldChange, y=-log10(padj_reduce)),
                                            color=restoplot$threshold)+
  geom_text_repel(data=right, aes(x=log2FoldChange, y=-log10(padj_reduce), label=GeneId),
                  segment.size = 0.5)+
  geom_text_repel(data=left, aes(x=log2FoldChange, y=-log10(padj_reduce), label=GeneId),
                  segment.size = 0.5)+
  geom_vline(aes(xintercept=1),color='brown')+
  geom_vline(aes(xintercept=-1),color='brown')+
  geom_hline(aes(yintercept=10), color='brown')+
  labs(title=paste0('Volcano plot: ', condition1, ' vs ' ,
                    condition2, '. UP: ',up.nrow,', DOWN: ',down.nrow), x='Effect size: log2(fold-change)',
       y='-log10(adjusted p-value)')+
  theme_bw()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12))+
  geom_text(label="Higher in DRG chronic pain", x=15, y=260, color='red', vjust=1)+
  geom_text(label="Higher in H9noc D28", x=-16, y=260, color='#007F00', vjust=1)+
  xlim(-20,20)+ylim(0,260) # may need to expand to -20,20 if the log2fold changes are that big
plot

# volcano plot, noci D28 versus DRGs control-----
library(data.table)
library(ggthemes)
library(ggrepel)
library(cowplot)

res <- fread('/Users/malleyce/Downloads/DRGs/DE/DE.H9noc_D28.vs.DRG_control.csv')
restoplot <- na.omit(res)
condition1 <- 'H9noc_D28'
condition2 <- 'DRG_control'

# set p to 1*10^-300 where padj=0 for plotting gene labels more easily.
restoplot[,padj_reduce:=padj]
restoplot[padj_reduce ==0,padj_reduce := (1*10^(-300))]

right <- subset(restoplot, ( ((log2FoldChange > 1) & (-log10(padj_reduce) > 10)))  ) #0
left <- subset(restoplot, ( ((log2FoldChange < -1) & (-log10(padj_reduce) > 10)))) #0
down.nrow <- nrow(left)
up.nrow <- nrow(right)

if(nrow(right)>15){
  right <- right[order(-log2FoldChange)]
  right <- right[1:15,]
}

if(nrow(left)>15){
  left <- left[order(log2FoldChange)]
  left <- left[1:15,]
}

restoplot[,threshold:=ifelse( ((-log10(padj_reduce) > 10) & log2FoldChange < -1),
                              '#007F00', 
                              ifelse( ((-log10(padj_reduce) > 10) & log2FoldChange > 1),
                                      'red', 'gray' ))]

plot <- ggplot(data=restoplot) + geom_point(aes(x=log2FoldChange, y=-log10(padj_reduce)),
                                            color=restoplot$threshold)+
  geom_text_repel(data=right, aes(x=log2FoldChange, y=-log10(padj_reduce), label=GeneId),
                  segment.size = 0.5)+
  geom_text_repel(data=left, aes(x=log2FoldChange, y=-log10(padj_reduce), label=GeneId),
                  segment.size = 0.5)+
  geom_vline(aes(xintercept=1),color='brown')+
  geom_vline(aes(xintercept=-1),color='brown')+
  geom_hline(aes(yintercept=10), color='brown')+
  labs(title=paste0('Volcano plot: ', condition1, ' vs ' ,
                    condition2, '. UP: ',up.nrow,', DOWN: ',down.nrow), x='Effect size: log2(fold-change)',
       y='-log10(adjusted p-value)')+
  theme_bw()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12))+
  geom_text(label="Higher in DRG control", x=15, y=133, color='red', vjust=1)+
  geom_text(label="Higher in H9noc D28", x=-16, y=133, color='#007F00', vjust=1)+
  xlim(-20,20)+ylim(0,133) # may need to expand to -20,20 if the log2fold changes are that big
plot

# H9nocD28 vs Axol, redo-----
res <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/DE/DE.H9noc_D28_new.vs.Axol_Noc.csv')
res
restoplot <- na.omit(res)
condition1 <- 'H9noc_D28'
condition2 <- 'Axol'
# set p to 1*10^-300 where padj=0 for plotting gene labels more easily.
restoplot[,padj_reduce:=padj]
restoplot[padj_reduce ==0,padj_reduce := (1*10^(-300))]

right <- subset(restoplot, ( ((log2FoldChange > 1) & (-log10(padj_reduce) > 10)))  ) #0
left <- subset(restoplot, ( ((log2FoldChange < -1) & (-log10(padj_reduce) > 10)))) #0
down.nrow <- nrow(left)
up.nrow <- nrow(right)

if(nrow(right)>15){
  right <- right[order(-log2FoldChange)]
  right <- right[1:15,]
}

if(nrow(left)>15){
  left <- left[order(log2FoldChange)]
  left <- left[1:15,]
}

restoplot[,threshold:=ifelse( ((-log10(padj_reduce) > 10) & log2FoldChange < -1),
                              '#007F00', 
                              ifelse( ((-log10(padj_reduce) > 10) & log2FoldChange > 1),
                                      'red', 'gray' ))]

plot <- ggplot(data=restoplot) + geom_point(aes(x=log2FoldChange, y=-log10(padj_reduce)),
                                            color=restoplot$threshold)+
  geom_text_repel(data=right, aes(x=log2FoldChange, y=-log10(padj_reduce), label=GeneId),
                  segment.size = 0.5)+
  geom_text_repel(data=left, aes(x=log2FoldChange, y=-log10(padj_reduce), label=GeneId),
                  segment.size = 0.5)+
  geom_vline(aes(xintercept=1),color='brown')+
  geom_vline(aes(xintercept=-1),color='brown')+
  geom_hline(aes(yintercept=10), color='brown')+
  labs(title=paste0('Volcano plot: ', condition1, ' vs ' ,
                    condition2, '. UP: ',up.nrow,', DOWN: ',down.nrow), x='Effect size: log2(fold-change)',
       y='-log10(adjusted p-value)')+
  theme_bw()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12))+
  geom_text(label="Higher in Axol", x=15, y=300, color='red', vjust=1)+
  geom_text(label="Higher in H9noc D28", x=-16, y=300, color='#007F00', vjust=1)+
  xlim(-20,20)+ylim(0,300) # may need to expand to -20,20 if the log2fold changes are that big
plot

# H9nocD56 vs Axol, set up DDS-----
#need to make a new dds.
# see /Users/malleyce/Downloads/DRGs/NociD0_D56_Axol/make_DDS_NociD0-D56_Axol.R

#load('/Users/malleyce/Downloads/DRGs/NociD0_D56_Axol/H9noc_D0-D56_Axol.DDS.RData')




#H9nocD56 vs Axol volcano----
res <- fread('/Users/malleyce/Downloads/DRGs/NociD0_D56_Axol/DE/DE.H9noc_D56.vs.Axol_Noc.csv')
res
restoplot <- na.omit(res)
condition1 <- 'H9noc_D56'
condition2 <- 'Axol'
# set p to 1*10^-300 where padj=0 for plotting gene labels more easily.
restoplot[,padj_reduce:=padj]
restoplot[padj_reduce ==0,padj_reduce := (1*10^(-300))]

right <- subset(restoplot, ( ((log2FoldChange > 1) & (-log10(padj_reduce) > 10)))  ) #0
left <- subset(restoplot, ( ((log2FoldChange < -1) & (-log10(padj_reduce) > 10)))) #0
down.nrow <- nrow(left)
up.nrow <- nrow(right)

if(nrow(right)>15){
  right <- right[order(-log2FoldChange)]
  right <- right[1:15,]
}

if(nrow(left)>15){
  left <- left[order(log2FoldChange)]
  left <- left[1:15,]
}

restoplot[,threshold:=ifelse( ((-log10(padj_reduce) > 10) & log2FoldChange < -1),
                              '#007F00', 
                              ifelse( ((-log10(padj_reduce) > 10) & log2FoldChange > 1),
                                      'red', 'gray' ))]

plot <- ggplot(data=restoplot) + geom_point(aes(x=log2FoldChange, y=-log10(padj_reduce)),
                                            color=restoplot$threshold)+
  geom_text_repel(data=right, aes(x=log2FoldChange, y=-log10(padj_reduce), label=GeneId),
                  segment.size = 0.5)+
  geom_text_repel(data=left, aes(x=log2FoldChange, y=-log10(padj_reduce), label=GeneId),
                  segment.size = 0.5)+
  geom_vline(aes(xintercept=1),color='brown')+
  geom_vline(aes(xintercept=-1),color='brown')+
  geom_hline(aes(yintercept=10), color='brown')+
  labs(title=paste0('Volcano plot: ', condition1, ' vs ' ,
                    condition2, '. UP: ',up.nrow,', DOWN: ',down.nrow), x='Effect size: log2(fold-change)',
       y='-log10(adjusted p-value)')+
  theme_bw()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12))+
  geom_text(label="Higher in Axol", x=17, y=300, color='red', vjust=1)+
  geom_text(label="Higher in H9noc D56", x=-16, y=300, color='#007F00', vjust=1)+
  xlim(-20,20)+ylim(0,300) # may need to expand to -20,20 if the log2fold changes are that big
plot
