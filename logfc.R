library(data.table)
library(DESeq2)
library(Hmisc)


norm <- as.data.frame(counts(dds, normalized=T))
norm$gene <- row.names(norm)
norm <- as.data.table(norm)
norm <- norm[,c(46, 1:45)]
norm

coldata <- as.data.frame(dds@colData)
write.csv(coldata, 'coldata.csv')


# auto increment sample replicates within the condition (treatment) group by _1, _2, etc.

sampleTable <- coldata
sampleTable$file <- row.names(sampleTable)
sampleTable <- as.data.table(sampleTable)
sampleTable[,replicate:= condition]
sampleTable[, id := seq_len(.N), by = replicate]
sampleTable[, id := rowid(replicate)]

sampleTable[,replicate:= paste0(replicate,'_',id)]
sampleTable[,id:=NULL]
sampleTable
sampleTable <- sampleTable[,c('file','condition','replicate')]


names(norm) <- c('gene',sampleTable$replicate)
names(norm)



###

norm[1:10,1:10]


norm.log <- norm

cols <- names(norm.log)[c(names(norm.log) %nin% c('gene','H9noc_D0_1','H9noc_D0_2','H9noc_D0_3'))]
norm.log[,AvgNorm_D0 := rowMeans(norm.log[,c('H9noc_D0_1','H9noc_D0_2','H9noc_D0_3')])]
norm.log[,AvgNorm_D0_plus1 := AvgNorm_D0 + 1]
norm.log[,(cols):= lapply(.SD, '+', 1), .SDcols=cols]
norm.log[,(cols):= lapply(.SD, '/', norm.log$AvgNorm_D0_plus1), .SDcols=cols]
norm.log[ , (cols) := lapply(.SD, log2), .SDcols = cols]

norm.log

names(norm.log)[2:4] <- c('H9noc_D0_1_exp', 'H9noc_D0_2_exp', 'H9noc_D0_3_exp')
names(norm.log)[5:46]
norm.log[,AvgNorm_D0:= NULL]
norm.log[,AvgNorm_D0_plus1:= NULL]
fwrite(norm.log, 'H9noc_DRGs_Log2FC.csv')

##



