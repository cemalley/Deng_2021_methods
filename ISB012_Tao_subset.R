library(RUVSeq) # http://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf
library(EDASeq)
library(RColorBrewer)
library(DESeq2)
library(data.table)
library(ComplexHeatmap)
library("scales")

setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/')
sampletable <- as.data.table(readxl::read_xlsx('sampletable.xlsx'))
sampletable

sampletable <- sampletable[-c(22),]

load('ISB012.Tao.DDS.pre_correction.RData')

counts <- as.data.frame(counts(dds, normalized=F))
counts <- counts[,-c(22)]
names(counts) <- sampletable$replicate
counts <- as.matrix(counts)

x <- as.factor(sampletable$condition)
set <- newSeqExpressionSet(counts,
                           phenoData = data.frame(x, row.names= c(sampletable$replicate) ))
set

colors <- brewer.pal(9, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)


spikes <- unique(c('ANAPC5', 'ANAPC15', 'ARID3B', 'ARL10', 'ATXN2', 'C16orf62', 'C3orf49', 'CCAR1', 'CCDC125', 'CCDC90B', 'CHFR', 'DHRSX', 'FRMD8', 'GGA1', 'HERC4', 'MKNK1', 'NASP', 'NME4', 'OTUB1', 'PMF1', 'POLR2B', 'POLR3A', 'POMK', 'PSMA3-AS1', 'PTPN14', 'RAPGEF6', 'REL', 'RRP1', 'RUNDC1', 'SAMD4B', 'SLC4A1AP', 'SLMAP', 'SMARCAL1', 'SNAP29', 'SNRNP200', 'SUPT4H1', 'TBC1D22A', 'THUMPD3-AS1', 'TSPOAP1-AS1', 'TUBGCP2', 'WDTC1', 'ZNF544','C1orf43','CHMP2A','EMC7','GPI','PSMB2','PSMB4','RAB7A','REEP5','SNRPD3','VCP','VPS29')) #https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-019-0538-z#Tab3

spikes <- spikes[spikes %in% row.names(counts)]

set2 <- RUVg(set, spikes, k=1)
pData(set2)

plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set2, col=colors[x], cex=1.2)


condition <- x

names(pData(set2)) <- c('condition', 'W_1')

dds <- DESeqDataSetFromMatrix(countData = counts(set2),
                              colData = pData(set2),
                              design = ~ W_1 + condition)
dds <- DESeq(dds)

save(dds, file='ISB012.Tao.DDS.after_correction_withoutH9D28New1.RData')

# september 4 set, with new d28 and with axol----

setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/')
sampletable <- as.data.table(readxl::read_xlsx('sampletable.xlsx', sheet=2))
sampletable

directory <- "/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/Countfiles_gene_symbols_intersection/"
dds <- DESeqDataSetFromHTSeqCount(sampletable,
                                  directory = directory,
                                  design = ~ condition)

dds <- DESeq(dds)


counts <- as.data.frame(counts(dds, normalized=F))
names(counts) <- sampletable$replicate
counts <- as.matrix(counts)

x <- as.factor(sampletable$condition)
set <- newSeqExpressionSet(counts,
                           phenoData = data.frame(x, row.names= c(sampletable$replicate) ))
set

colors <- brewer.pal(9, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)


spikes <- unique(c('ANAPC5', 'ANAPC15', 'ARID3B', 'ARL10', 'ATXN2', 'C16orf62', 'C3orf49', 'CCAR1', 'CCDC125', 'CCDC90B', 'CHFR', 'DHRSX', 'FRMD8', 'GGA1', 'HERC4', 'MKNK1', 'NASP', 'NME4', 'OTUB1', 'PMF1', 'POLR2B', 'POLR3A', 'POMK', 'PSMA3-AS1', 'PTPN14', 'RAPGEF6', 'REL', 'RRP1', 'RUNDC1', 'SAMD4B', 'SLC4A1AP', 'SLMAP', 'SMARCAL1', 'SNAP29', 'SNRNP200', 'SUPT4H1', 'TBC1D22A', 'THUMPD3-AS1', 'TSPOAP1-AS1', 'TUBGCP2', 'WDTC1', 'ZNF544','C1orf43','CHMP2A','EMC7','GPI','PSMB2','PSMB4','RAB7A','REEP5','SNRPD3','VCP','VPS29')) #https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-019-0538-z#Tab3

spikes <- spikes[spikes %in% row.names(counts)]

set2 <- RUVg(set, spikes, k=1)
pData(set2)

plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set2, col=colors[x], cex=1.2)


condition <- x

names(pData(set2)) <- c('condition', 'W_1')

dds <- DESeqDataSetFromMatrix(countData = counts(set2),
                              colData = pData(set2),
                              design = ~ W_1 + condition)
dds <- DESeq(dds)

save(dds, file='ISB012.Tao.DDS.after_correction_withH9D28New_Axol.RData')

# heatmap----
#genelist <- c('POU5F1', 'NANOG', 'SOX2', 'SOX21', 'OTX2', 'LHX2', 'PAX3', 'NEUROG1', 'TFAP2A', 'TFAP2B', 'SOX10', 'POU4F1', 'ISL1', 'DRGX', 'SLC17A6', 'PRPH', 'RUNX1', 'TAC1', 'RET', 'P2RX3', 'TRPV1', 'SCN9A', 'SCN10A', 'SCN11A', 'NTRK2', 'NRG2', 'PDGFRA', 'BDNF', 'NGF', 'SOX17', 'GATA4', 'GATA6', 'TBXT', 'EOMES')

genelist <- c('POU5F1', 'NANOG', 'SOX2', 'SOX21', 'OTX2', 'PAX3', 'NEUROG1', 'TFAP2A', 'TFAP2B', 'SOX10', 'POU4F1', 'ISL1', 'DRGX', 'SLC17A6', 'PRPH', 'RUNX1', 'TAC1', 'RET', 'P2RX3', 'TRPV1', 'SCN9A', 'SCN10A', 'SCN11A', 'NTRK2', 'NRG2', 'BDNF', 'NGF')


dds.stabilized <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")

dds.genes <- row.names(dds.stabilized)
found <- dds.genes[dds.genes %in% genelist]
length(found)
length(genelist)
genelist[!genelist %in% found] #

mat <- subset(assay(dds.stabilized), row.names(dds.stabilized) %in% genelist)

coldata <- as.data.frame(dds@colData)

mat <- as.data.frame(mat)
mat

# rearrange genes.-------
#mat <- as.data.frame(readxl::read_xlsx('mat_for_heatmap_subset.xlsx'))
#row.names(mat) <- mat$GeneId
#mat <- mat[-1]
#mat
#names(mat)<- gsub('_',' ',coldata$condition)

mat.dt <- data.table(GeneId = genelist, order=c(1:length(genelist)))
mat.dt

mat.original <- as.data.frame(mat)
mat.original$GeneId <- row.names(mat.original)
mat.original <- as.data.table(mat.original)

mat.dt <- merge(mat.original,mat.dt, by='GeneId')
mat.dt <- mat.dt[order(order)]

mat.dt

mat <- as.data.frame(mat.dt)
row.names(mat) <- mat$GeneId
mat <- mat[-23]
mat <- mat[-1]

cols.use <- colorRampPalette(colors=rev(brewer.pal(11,"RdBu")))(100) # reversed RdBu, creates Blue-white-red

scaled_mat <- t(scale(t(mat)))

rescaled_mat <- rescale(scaled_mat, to=c(-2,2))

Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = T, show_row_names = T,
        cluster_rows = FALSE, cluster_columns = FALSE,
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))


# january 8 2020 set, old day 0-21, plus new d28----
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/')
sampletable <- as.data.table(readxl::read_xlsx('sampletable.xlsx', sheet=3))
sampletable

directory <- "/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/Countfiles_gene_symbols_intersection/"
dds <- DESeqDataSetFromHTSeqCount(sampletable,
                                  directory = directory,
                                  design = ~ condition)

dds <- DESeq(dds)


counts <- as.data.frame(counts(dds, normalized=F))
names(counts) <- sampletable$replicate
counts <- as.matrix(counts)

x <- as.factor(sampletable$condition)
set <- newSeqExpressionSet(counts,
                           phenoData = data.frame(x, row.names= c(sampletable$replicate) ))
set

colors <- brewer.pal(9, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)


spikes <- unique(c('ANAPC5', 'ANAPC15', 'ARID3B', 'ARL10', 'ATXN2', 'C16orf62', 'C3orf49', 'CCAR1', 'CCDC125', 'CCDC90B', 'CHFR', 'DHRSX', 'FRMD8', 'GGA1', 'HERC4', 'MKNK1', 'NASP', 'NME4', 'OTUB1', 'PMF1', 'POLR2B', 'POLR3A', 'POMK', 'PSMA3-AS1', 'PTPN14', 'RAPGEF6', 'REL', 'RRP1', 'RUNDC1', 'SAMD4B', 'SLC4A1AP', 'SLMAP', 'SMARCAL1', 'SNAP29', 'SNRNP200', 'SUPT4H1', 'TBC1D22A', 'THUMPD3-AS1', 'TSPOAP1-AS1', 'TUBGCP2', 'WDTC1', 'ZNF544','C1orf43','CHMP2A','EMC7','GPI','PSMB2','PSMB4','RAB7A','REEP5','SNRPD3','VCP','VPS29')) #https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-019-0538-z#Tab3

spikes <- spikes[spikes %in% row.names(counts)]

set2 <- RUVg(set, spikes, k=1)
pData(set2)

plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set2, col=colors[x], cex=1.2)


condition <- x

names(pData(set2)) <- c('condition', 'W_1')

dds <- DESeqDataSetFromMatrix(countData = counts(set2),
                              colData = pData(set2),
                              design = ~ W_1 + condition)
dds <- DESeq(dds)

save(dds, file='ISB012.Tao.DDS.after_correction.RData')

# heatmap 1/24/20-----
load('ISB012.Tao.DDS.after_correction.RData')

genelist <- fread('markers_to_add.txt', header=F, nrows = 56)
genelist

genelist.previous <- c('POU5F1', 'NANOG', 'SOX2', 'SOX21', 'OTX2', 'PAX3', 'NEUROG1', 'TFAP2A', 'TFAP2B', 'SOX10', 'POU4F1', 'ISL1', 'DRGX', 'SLC17A6', 'PRPH', 'RUNX1', 'TAC1', 'RET', 'P2RX3', 'TRPV1', 'SCN9A', 'SCN10A', 'SCN11A', 'NTRK2', 'NRG2', 'BDNF', 'NGF')

genelist.merged <- unique(c(genelist.previous, genelist$V1))

genelist.merged

dds.stabilized <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")

dds.genes <- row.names(dds.stabilized)
found <- dds.genes[dds.genes %in% genelist.merged]
length(found)
length(genelist.merged)
genelist.merged[!genelist.merged %in% found] #

mat <- subset(assay(dds.stabilized), row.names(dds.stabilized) %in% genelist.merged)

coldata <- as.data.frame(dds@colData)

mat <- as.data.frame(mat)
mat

cols.use <- colorRampPalette(colors=rev(brewer.pal(11,"RdBu")))(100) # reversed RdBu, creates Blue-white-red

scaled_mat <- t(scale(t(mat)))

rescaled_mat <- rescale(scaled_mat, to=c(-2,2))

Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = T, show_row_names = T,
        cluster_rows = TRUE, cluster_columns = FALSE,
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))

Heatmap(as.matrix(scaled_mat), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = T, show_row_names = T,
        cluster_rows = TRUE, cluster_columns = FALSE,
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))

Heatmap(as.matrix(mat), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = T, show_row_names = T,
        cluster_rows = TRUE, cluster_columns = FALSE,
        heatmap_legend_param = list(legend_height = unit(4, "cm"), title='Normalized\nexpression'))

raw.counts <- as.data.frame(counts(dds, normalized=F))
names(raw.counts)
raw.counts$GeneId <- row.names(raw.counts)
raw.counts <- as.data.table(raw.counts)
raw.counts <- subset(raw.counts, select=c(19, 1:18))
raw.counts


# trying column wise scaling

col_scaled_mat <- scale(mat)
col_rescaled_mat <- rescale(col_scaled_mat, to=c(-2,2))
Heatmap(as.matrix(col_rescaled_mat), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = T,
        show_row_names = T,
        cluster_rows = TRUE, cluster_columns = FALSE,
        heatmap_legend_param = list(legend_height = unit(4, "cm"),
                                    title='Scaled\nexpression'))


# heatmap 1/27/20----
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/')
genelist <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/Marker genes in category.xlsx'))
genelist

dds.stabilized <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")

dds.genes <- row.names(dds.stabilized)
found <- dds.genes[dds.genes %in% genelist$GeneID]
length(found)
length(genelist$GeneID)
genelist$GeneID[!genelist$GeneID %in% found] #

mat <- subset(assay(dds.stabilized), row.names(dds.stabilized) %in% genelist$GeneID)

coldata <- as.data.frame(dds@colData)

mat <- as.data.frame(mat)
mat

cols.use <- colorRampPalette(colors=rev(brewer.pal(11,"RdBu")))(100) # reversed RdBu, creates Blue-white-red

scaled_mat <- t(scale(t(mat)))
rescaled_mat <- rescale(scaled_mat, to=c(-2,2))

Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = T,
        show_row_names = T,
        cluster_rows = TRUE, cluster_columns = FALSE,
        heatmap_legend_param = list(legend_height = unit(4, "cm"),
                                    title='Scaled\nexpression'))

col_scaled_mat <- scale(mat)
col_rescaled_mat <- rescale(col_scaled_mat, to=c(-2,2))
Heatmap(as.matrix(col_rescaled_mat), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = T,
        show_row_names = T,
        cluster_rows = TRUE, cluster_columns = FALSE,
        heatmap_legend_param = list(legend_height = unit(4, "cm"),
                                    title='Scaled\nexpression'))

# reorder
mat.dt <- data.table(GeneId = genelist$GeneID, order=c(1:length(genelist$GeneID)))
mat.dt

mat.original <- as.data.frame(mat)
mat.original$GeneId <- row.names(mat.original)
mat.original <- as.data.table(mat.original)

mat.dt <- merge(mat.original,mat.dt, by='GeneId')
mat.dt <- mat.dt[order(order)]

mat.dt

mat <- as.data.frame(mat.dt)
row.names(mat) <- mat$GeneId
mat <- mat[-23]
mat <- mat[-1]

col_scaled_mat <- scale(mat)
col_rescaled_mat <- rescale(col_scaled_mat, to=c(-2,2))

row_scaled_mat <- t(scale(t(mat)))
row_rescaled_mat <- rescale(row_scaled_mat, to=c(-2,2))

Heatmap(as.matrix(col_rescaled_mat), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = T,
        show_row_names = T,
        cluster_rows = TRUE, cluster_columns = FALSE,
        heatmap_legend_param = list(legend_height = unit(4, "cm"),
                                    title='Scaled\nexpression'),
        split=c(rep('Pluripotency', 2),
                rep('Neural ectoderm', 3),
                rep('Neural crest', 5),
                rep('Nociceptor', 13),
                rep('Peptidergic', 3),
                rep('Nonpeptidergic', 4),
                rep('Neurofilament', 2)), gap=unit(3, "mm")  )


Heatmap(as.matrix(mat), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = T,
        show_row_names = T,
        cluster_rows = TRUE, cluster_columns = FALSE,
        heatmap_legend_param = list(legend_height = unit(4, "cm"),
                                    title='Normalized\nexpression'),
        split=c(rep('Pluripotency', 2),
                rep('Neural ectoderm', 3),
                rep('Neural crest', 5),
                rep('Nociceptor', 13),
                rep('Peptidergic', 3),
                rep('Nonpeptidergic', 4),
                rep('Neurofilament', 2)), gap=unit(3, "mm"))


Heatmap(as.matrix(row_rescaled_mat), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = T,
        show_row_names = T,
        cluster_rows = FALSE, cluster_row_slices=FALSE, cluster_columns = FALSE,
        heatmap_legend_param = list(legend_height = unit(4, "cm"),
                                    title='Row Z-score'),
        split= factor(c(rep('Pluripotency', 2),
                rep('Neural ectoderm', 3),
                rep('Neural crest', 5),
                rep('Nociceptor', 13),
                rep('Peptidergic', 3),
                rep('Nonpeptidergic', 3),
                rep('Neurofilament', 2)), levels=c('Pluripotency',
                                                   'Neural ectoderm',
                                                   'Neural crest',
                                                   'Nociceptor',
                                                   'Peptidergic',
                                                   'Nonpeptidergic',
                                                   'Neurofilament')), gap=unit(2, "mm"))

# february 13, 2020----
load("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/ISB012.Tao.DDS.after_correction_withoutH9D28New1.RData")

View(as.data.frame(dds@colData))
raw.counts <- counts(dds, normalized=F)
raw.counts.nociD8 <- subset(raw.counts, select=c('H9noc_D0_1','H9noc_D0_1','H9noc_D0_1','H9noc_D8_1', 'H9noc_D8_2','H9noc_D8_3'))
raw.counts.nociD8 <- as.data.frame(raw.counts.nociD8)
raw.counts.nociD8
fwrite(raw.counts.nociD8, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/ISB012_TaoNociD8_rawcounts.csv', row.names = T, col.names = T, quote=F, sep='\t')

norm.counts <-  counts(dds, normalized=T)
norm.counts.nociD8 <- subset(norm.counts, select=c('H9noc_D0_1','H9noc_D0_1','H9noc_D0_1','H9noc_D8_1', 'H9noc_D8_2','H9noc_D8_3'))
norm.counts.nociD8 <- as.data.frame(norm.counts.nociD8)
norm.counts.nociD8
fwrite(norm.counts.nociD8, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/ISB012_TaoNociD8_normcounts.csv', row.names = T, col.names = T, quote=F, sep='\t')

# volcano plot for d28 new vs axol, may 12 2020----
library(data.table)
library(DESeq2)
library(scales)
library(ggthemes)
library(ggplot2)
res <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/DE/DE.H9noc_D28_new.vs.Axol_Noc.csv')

restoplot <- na.omit(res)
restoplot[,padj := as.numeric(padj)]

restoplot[,threshold:=ifelse(padj <=(1*(10^(-2))) & log2FoldChange > 2, '#007F00', 
                             ifelse( padj <=(1*(10^(-2))) & log2FoldChange < -2, 'red', 'gray' ))]

restoplot[,padj_reduce := ifelse(padj <= (1*(10^(-10))), (1*(10^(-10))), padj)]

right <- subset(restoplot, threshold=='#007F00' & log2FoldChange > 2)
left <- subset(restoplot, threshold=='red' & log2FoldChange < 2)

condition1 <- 'H9noc_D28_new'
condition2 <- 'Axol_Noc'

# no gene labels:
attach(restoplot)
class(restoplot$threshold)
factor(restoplot$threshold)

restoplot[,threshold:= factor(threshold, levels=c('gray','red', '#007F00'))]

restoplot[,log2FoldChange_adj:= ifelse(log2FoldChange >=25, 25, log2FoldChange)]
restoplot[,log2FoldChange_adj:= ifelse(log2FoldChange_adj <= -25, -25, log2FoldChange_adj)]

ggplot(data=restoplot) + geom_point(aes(x=log2FoldChange_adj,
    y=-log10(padj_reduce)), size=2,color=restoplot$threshold)+
  labs(title=paste0('Volcano plot: ', condition1, ' vs ' ,
       condition2), x='Effect size: log2(fold-change)',
       y='-log10(adjusted p-value)')+
   theme_hc()+ scale_x_continuous(breaks=c(-25, -15,0, 15, 25)) +
        scale_y_continuous(breaks=c(0,2,4,6,8,10))+
        theme(axis.title = element_text(size=15),
              axis.text = element_text(size=12))

restoplot[threshold =='red',] #5 down
restoplot[threshold =='#007F00',] #56 up

# with gene labels:
res <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/DE/DE.H9noc_D28_new.vs.Axol_Noc.csv')

restoplot <- na.omit(res)
condition1 <- 'H9noc_D28'
condition2 <- 'Axol_noc'

# set p to 1*10^-300 where padj=0 for plotting gene labels more easily.
restoplot[,padj_reduce:=padj]
restoplot[padj_reduce ==0,padj_reduce := (1*10^(-300))]

right <- subset(restoplot, ((log2FoldChange > 1) & (-log10(padj_reduce) > 50))) #236
left <- subset(restoplot, ( ((log2FoldChange < -1) & (-log10(padj_reduce) > 50)))) #210

if(nrow(right)>15){
        right <- right[order(-log2FoldChange)]
        right <- right[1:15,]
}

if(nrow(left)>15){
        left <- left[order(log2FoldChange)]
        left <- left[1:15,]
}

restoplot[,threshold:=ifelse((-log10(padj_reduce) > 50) & log2FoldChange > 1, '#007F00', 
                             ifelse( (-log10(padj_reduce) > 50) & log2FoldChange < -1, 'red', 'gray' ))]

plot <- ggplot(data=restoplot) + geom_point(aes(x=log2FoldChange, y=-log10(padj_reduce)),
                                            color=restoplot$threshold)+
        geom_text_repel(data=right, aes(x=log2FoldChange, y=-log10(padj_reduce), label=GeneId),
                        segment.size = 0.5)+
        geom_text_repel(data=left, aes(x=log2FoldChange, y=-log10(padj_reduce), label=GeneId),
                        segment.size = 0.5)+
        geom_vline(aes(xintercept=1),color='brown')+
        geom_vline(aes(xintercept=-1),color='brown')+
        geom_hline(aes(yintercept=50), color='brown')+
        labs(title=paste0('Volcano plot: ', condition1, ' vs ' ,
                          condition2, '. UP: 236, DOWN: 210'), x='Effect size: log2(fold-change)',
             y='-log10(adjusted p-value)')+
        theme_bw()+
        theme(axis.title = element_text(size=15),
              axis.text = element_text(size=12))+
        xlim(-15,15)+ylim(0,320) # may need to expand to -20,20 if the log2fold changes are that big
plot


# enrichr for d28 vs axol-----
res <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/DE/DE.H9noc_D28_new.vs.Axol_Noc.csv')
res <- res[padj <= 0.001,]
res <- res[order(-log2FoldChange)]
res
up_axol <- res[log2FoldChange>0,GeneId]
cat(up_axol[1:100], sep='\n')

res <- res[order(log2FoldChange)]
up_h9nocd28 <- res[log2FoldChange<0,GeneId]
cat(up_h9nocd28[1:100], sep='\n')

# heatmap of ion channels-----
load('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/ISB012.Tao.DDS.after_correction_withH9D28New_Axol.RData')
library(DESeq2)
library(ComplexHeatmap)
library(data.table)
library(RColorBrewer)
library(scales)

dds.subset <- dds[,colData(dds)$condition %in% c('H9noc_D28_new','Axol_Noc')]

dds.stabilized <- varianceStabilizingTransformation(dds.subset, blind = TRUE, fitType = "parametric")

genelist <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/marker_sets/Ion channels group_updated.xlsx'))
#genelist <- genelist[order(`Approved symbol`)]
#genelist <- genelist[!duplicated(`Approved symbol`)]
unique(genelist$`Group name`)
# [1] "Chloride channels"          "Porins"                    
# [3] "Sodium channels"            "Calcium channels"          
# [5] "Voltage-gated ion channels" "Ligand gated ion channels" 
# [7] "Gap junction proteins"      "Potassium channels" 

chloride.channel <- genelist[`Group name`=='Chloride channels',c(`Approved symbol`)]
porins <- genelist[`Group name`=='Porins',c(`Approved symbol`)]
sodium.channel <- genelist[`Group name`=='Sodium channels',c(`Approved symbol`)]
calcium.channel <- genelist[`Group name`=='Calcium channels',c(`Approved symbol`)]
voltage <- genelist[`Group name`=='Voltage-gated ion channels',c(`Approved symbol`)]
ligand <- genelist[`Group name`=='Ligand gated ion channels',c(`Approved symbol`)]
gap <- genelist[`Group name`=='Gap junction proteins',c(`Approved symbol`)]
potassium <- genelist[`Group name`=='Potassium channels',c(`Approved symbol`)]

#chloride channel
genes <- unique(potassium)
dds.genes <- row.names(dds.stabilized)
found <- dds.genes[dds.genes %in% genes]
length(found)
length(genes)
genes[!genes %in% found] #[1] "GABRQ"  "GJA5"   "KCNJ18"
mat <- subset(assay(dds.stabilized), row.names(dds.stabilized) %in% genes)
coldata <- as.data.frame(dds@colData)
mat <- as.data.frame(mat)
names(mat)<- c(rep('H9noc_D28', 3), rep('Axol_Noc',3))
mat
cols.use <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100) # reversed RdBu, creates Blue-white-red
scaled_mat <- t(scale(t(mat)))
rescaled_mat <- rescale(scaled_mat, to=c(-2,2))
rescaled_mat  <- na.omit(rescaled_mat)
# h1 <- Heatmap(as.matrix(rescaled_mat),column_title = "Potassium channels", row_names_side = "right",
#               column_names_side = "top",col = cols.use, show_column_names = T,
#               cluster_rows = TRUE, cluster_columns = FALSE,
#               heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))
h2 <- Heatmap(as.matrix(mat),column_title = "Potassium channels", row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = T,
              cluster_rows = TRUE, cluster_columns = FALSE,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Normalized expression'))

#chloride.heatmap <- h1
#porins.heatmap <- h1
#sodium.channel.ht <- h1
#calcium.channel.ht <- h1
#voltage.ht <- h1
#ligand.ht <- h1
#gap.ht <- h1
#potassium.ht <- h1

chloride.normexp.ht <- h2
porins.normexp.ht <- h2
sodium.normexp.ht <- h2
calcium.normexp.ht <- h2
voltage.normexp.ht <- h2
ligand.normexp.ht <- h2
gap.normexp.ht <- h2
potassium.normexp.ht <- h2




draw(chloride.heatmap)
draw(porins.heatmap)
draw(sodium.channel.ht)
draw(calcium.channel.ht)
draw(voltage.ht)
draw(ligand.ht)
draw(gap.ht)
draw(potassium.ht)


chloride.normexp.ht
porins.normexp.ht
sodium.normexp.ht
calcium.normexp.ht
voltage.normexp.ht
ligand.normexp.ht
gap.normexp.ht
potassium.normexp.ht


# venn diagram for expression past a cutoff----
# let's set cutoff for "expressed" to 10 for norm counts.

genes.present <- unique(c(chloride.channel, porins, sodium.channel,
                          calcium.channel, voltage, ligand,
                          gap, potassium))[unique(c(chloride.channel, porins, sodium.channel,
                        calcium.channel, voltage, ligand,
                        gap, potassium)) %in% row.names(counts(dds.subset))]

genelist.dt <- data.table(GeneId=c(chloride.channel, porins, sodium.channel,
                                   calcium.channel, voltage, ligand,
                                   gap, potassium),
                          Category=c(rep('Chloride channels', length(chloride.channel)),
                                     rep('Porins', length(porins)),
                                     rep('Sodium channels', length(sodium.channel)),
                                     rep('Calcium channels', length(calcium.channel)),
                                     rep('Voltage-gated ion channels', length(voltage)),
                                     rep('Ligand-gated ion channels', length(ligand)),
                                     rep('Gap junciton proteins', length(gap)),
                                     rep('Potassium channels', length(potassium))))

genelist.dt <- subset(genelist.dt, genelist.dt$GeneId %in% c(row.names(counts(dds.subset))))
genelist.dt #441 genes with duplciates.

dds.subset
dds.subset.ions <- subset(dds.subset, row.names(counts(dds.subset)) %in% genelist.dt$GeneId) #325 unique genes.

normexp <- as.data.frame(counts(dds.subset.ions, normalized=T))
normexp[1:5,]
normexp$GeneId <- row.names(normexp)

nrow(merge(genelist.dt, normexp, by='GeneId'))

genelist.dt <- merge(genelist.dt, normexp, by='GeneId')
genelist.dt[,H9noc_avg := rowMeans(.SD),.SDcols=c('H9noc_D28_new_1','H9noc_D28_new_2','H9noc_D28_new_3')]
genelist.dt[,Axol_avg := rowMeans(.SD),.SDcols=c('Axol_Noc_1','Axol_Noc_2','Axol_Noc_3')]
genelist.dt

genelist.dt[H9noc_avg >=10 & Category=='Porins',.N, ]

# get category sizes:
#dat[, .(count = .N, var = sum(VAR)), by = MNTH]
genelist.dt[, .(count = .N), by = Category]

# overall venn for all ion channels

unique.ions <- genelist.dt[!duplicated(GeneId)]
unique.ions

nrow(unique.ions[H9noc_avg > 10,]) #235/325
nrow(unique.ions[Axol_avg > 10,]) #193/325


cat(unique.ions[H9noc_avg > 10,c('GeneId')]$GeneId, sep='\n')

cat(unique.ions[Axol_avg > 10,c('GeneId')]$GeneId, sep='\n')
#

# scatterplot of all genes-----
str(normexp)
normexp[,H9noc_avg := rowMeans(.SD),.SDcols=c('H9noc_D28_new_1','H9noc_D28_new_2','H9noc_D28_new_3')]
normexp[,Axol_avg := rowMeans(.SD),.SDcols=c('Axol_Noc_1','Axol_Noc_2','Axol_Noc_3')]


cowplot::plot_grid(ggplot(data=normexp, aes(x=H9noc_avg, y=Axol_avg))+geom_point()+
  geom_abline(slope=1, intercept = 0)+theme_bw()+
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=15))+
    ylim(0,220000)+xlim(0,220000)+
    labs(x='H9noc D28', y='Axol noc', title='Average normalized expression in SCTL and Axol nociceptors'))

normexp.label <- normexp[H9noc_avg>60000 | Axol_avg>60000,]

cowplot::plot_grid(ggplot(data=normexp, aes(x=H9noc_avg, y=Axol_avg))+geom_point()+
                     geom_abline(slope=1, intercept = 0)+theme_bw()+
                     theme(axis.text=element_text(size=15),
                           axis.title=element_text(size=15))+
                     ylim(0,220000)+xlim(0,220000)+
                     labs(x='H9noc D28', y='Axol noc', title='Average normalized expression in SCTL and Axol nociceptors')+
                     geom_text_repel(data=normexp.label,
                            aes(x=H9noc_avg, y=Axol_avg, label=GeneId))
                     )

#

# full transcriptome venn----
nrow(normexp[H9noc_avg > 10,]) #15301/32180
nrow(normexp[Axol_avg > 10,]) #14224/32180


expressed_H9noc <- normexp[H9noc_avg > 10,c('GeneId')]

expressed_Axol <- normexp[Axol_avg > 10,c('GeneId')]

fwrite(expressed_H9noc,file='/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/Venn/Ion_channels_over_10/expressed_h9noc.csv', col.names=F)
fwrite(expressed_Axol,file='/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/Venn/Ion_channels_over_10/expressed_axol.csv', col.names=F)

#
