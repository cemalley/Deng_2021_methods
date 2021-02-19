library(reshape2)
library(data.table)
library(DESeq2)
library(ComplexHeatmap)
library(RColorBrewer)

setwd('/Users/malleyce/Downloads/DRGs/')
genelist <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/Marker genes in category.xlsx'))
genelist

load('DRGs.Nociceptors.DDS.RData')

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
mat <- mat[-47]
mat <- mat[-1]

col_scaled_mat <- scale(mat)
col_rescaled_mat <- rescale(col_scaled_mat, to=c(-2,2))

row_scaled_mat <- t(scale(t(mat)))
row_rescaled_mat <- rescale(row_scaled_mat, to=c(-2,2))

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
                                                           'Neurofilament')), gap=unit(1, "mm"))


# heatmap of ion channels-----
load('DRGs.Nociceptors.DDS.RData')
library(DESeq2)
library(ComplexHeatmap)
library(data.table)
library(RColorBrewer)
library(scales)

dds.subset <- dds[,colData(dds)$condition %in% c('H9noc_D28','H9noc_D56','DRG_control','DRG_chronic_pain')]

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
genes[!genes %in% found] #[1]  "KCNJ18"
mat <- subset(assay(dds.stabilized), row.names(dds.stabilized) %in% genes)
coldata <- as.data.frame(dds.subset@colData)
mat <- as.data.frame(mat)
names(mat)<- c(rep('H9noc_D28', 3), rep('H9noc_D56',3),
               rep('DRG_control',3), rep('DRG_chronic_pain',21))
mat
cols.use <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100) # reversed RdBu, creates Blue-white-red
# scaled_mat <- t(scale(t(mat)))
# rescaled_mat <- rescale(scaled_mat, to=c(-2,2))
# rescaled_mat  <- na.omit(rescaled_mat)
# h1 <- Heatmap(as.matrix(rescaled_mat),column_title = "Potassium channels", row_names_side = "right",
#               column_names_side = "top",col = cols.use, show_column_names = T,
#               cluster_rows = TRUE, cluster_columns = FALSE,
#               heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))
h2 <- Heatmap(as.matrix(mat),column_title = "Potassium channels", row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = T,
              cluster_rows = TRUE, cluster_columns = FALSE,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Normalized expression'))


chloride.normexp.ht <- h2
porins.normexp.ht <- h2
sodium.normexp.ht <- h2
calcium.normexp.ht <- h2
voltage.normexp.ht <- h2
ligand.normexp.ht <- h2
gap.normexp.ht <- h2
potassium.normexp.ht <- h2


chloride.normexp.ht
porins.normexp.ht
sodium.normexp.ht
calcium.normexp.ht
voltage.normexp.ht
ligand.normexp.ht
gap.normexp.ht
potassium.normexp.ht


