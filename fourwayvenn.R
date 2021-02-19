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
genelist.dt #439 genes with duplciates. 323 unique.

dds.subset
dds.subset.ions <- subset(dds.subset, row.names(counts(dds.subset)) %in% genelist.dt$GeneId) #323 unique genes.

normexp <- as.data.frame(counts(dds.subset.ions, normalized=T))
normexp[1:5,]
normexp$GeneId <- row.names(normexp)

nrow(merge(genelist.dt, normexp, by='GeneId'))

genelist.dt <- merge(genelist.dt, normexp, by='GeneId')
names(genelist.dt) <- c('GeneId', 'Category', 'H9noc_D28_1', 'H9noc_D28_2', 'H9noc_D28_3', 'H9noc_D56_1', 'H9noc_D56_2', 'H9noc_D56_3', 'DRG_control_1', 'DRG_control_2', 'DRG_control_3', 'DRG_chronic_pain_1', 'DRG_chronic_pain_2', 'DRG_chronic_pain_3', 'DRG_chronic_pain_4', 'DRG_chronic_pain_5', 'DRG_chronic_pain_6', 'DRG_chronic_pain_7', 'DRG_chronic_pain_8', 'DRG_chronic_pain_9', 'DRG_chronic_pain_10', 'DRG_chronic_pain_11', 'DRG_chronic_pain_12', 'DRG_chronic_pain_13', 'DRG_chronic_pain_14', 'DRG_chronic_pain_15', 'DRG_chronic_pain_16', 'DRG_chronic_pain_17', 'DRG_chronic_pain_18', 'DRG_chronic_pain_19', 'DRG_chronic_pain_20', 'DRG_chronic_pain_21')
genelist.dt[,H9nocD28_avg := rowMeans(.SD),.SDcols=c('H9noc_D28_1','H9noc_D28_2','H9noc_D28_3')]
genelist.dt[,H9nocD56_avg := rowMeans(.SD),.SDcols=c('H9noc_D56_1','H9noc_D56_2','H9noc_D56_3')]
genelist.dt[,DRGcontrol_avg := rowMeans(.SD),.SDcols=c('DRG_control_1','DRG_control_2','DRG_control_3')]
genelist.dt[,DRGpain_avg := rowMeans(.SD),.SDcols=c('DRG_chronic_pain_1', 'DRG_chronic_pain_2', 'DRG_chronic_pain_3', 'DRG_chronic_pain_4', 'DRG_chronic_pain_5', 'DRG_chronic_pain_6', 'DRG_chronic_pain_7', 'DRG_chronic_pain_8', 'DRG_chronic_pain_9', 'DRG_chronic_pain_10', 'DRG_chronic_pain_11', 'DRG_chronic_pain_12', 'DRG_chronic_pain_13', 'DRG_chronic_pain_14', 'DRG_chronic_pain_15', 'DRG_chronic_pain_16', 'DRG_chronic_pain_17', 'DRG_chronic_pain_18', 'DRG_chronic_pain_19', 'DRG_chronic_pain_20', 'DRG_chronic_pain_21')]

genelist.dt[H9nocD28_avg >=10 & Category=='Porins',.N, ]

# get category sizes:
#dat[, .(count = .N, var = sum(VAR)), by = MNTH]
genelist.dt[, .(count = .N), by = Category]

# overall venn for all ion channels

unique.ions <- genelist.dt[!duplicated(GeneId)]
unique.ions

# nrow(unique.ions[H9noc_avg > 10,]) #235/325
# nrow(unique.ions[Axol_avg > 10,]) #193/325


H9nocD28_genes <- unique.ions[H9nocD28_avg > 10,c('GeneId')]$GeneId

H9nocD56_genes <- unique.ions[H9nocD56_avg > 10,c('GeneId')]$GeneId

DRGcontrol_genes <- unique.ions[DRGcontrol_avg > 10,c('GeneId')]$GeneId

DRGpain_genes <- unique.ions[DRGpain_avg > 10,c('GeneId')]$GeneId


A <- H9nocD28_genes
B <- H9nocD56_genes
C <- DRGcontrol_genes
D <- DRGpain_genes

length(A %in% B) #235


# four way venn-----
library(gplots)
library(gplots)
test1 <- c("dog", "cat", "monkey", "fish", "cow", "frog")
test2 <- c("cat", "frog", "aardvark", "monkey", "cow", "lizard", "bison", "goat")
test3 <- c("whale", "cat", "cow", "dog", "worm") 
test4 <- c("dog", "bird", "plant", "fly", "cow", "horse", "goat")

venn(list(A=test1,B=test2,C=test3,D=test4))

plot <- venn(list(A=H9nocD28_genes, B=H9nocD56_genes,C=DRGcontrol_genes,D=DRGpain_genes))
svg(plot, filename = "/Users/malleyce/Downloads/DRGs/NociD0_D56_Axol/Rplot.svg")
ggsave(filename = "testing.svg")
pdf(plot, 'venn.pdf')
dev.off()
venn(list(A=H9nocD28_genes, B=DRGpain_genes))
venn(list(A=H9nocD28_genes, B=DRGcontrol_genes))
venn(list(A=H9nocD56_genes, B=DRGpain_genes))
venn(list(A=H9nocD56_genes, B=DRGcontrol_genes))


xx.1 <- list(A = H9nocD28_genes, 
             B = H9nocD56_genes, 
             C = DRGcontrol_genes, 
             D = DRGpain_genes)

library(VennDiagram) # only working on biowulf right now.
overlap <- calculate.overlap(xx.1)
names(overlap) <- c("a1234", "a123", "a124", "a134", "a234", "a12", "a13", "a14", "a23", "a24", "a34", "a1", "a2", "a3", "a4") #correctly renames the sections