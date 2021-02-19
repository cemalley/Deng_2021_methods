library(data.table)
library(DESeq2)
library(readxl)
library(stringr)

samples_DRG_paired <- paste0('SRR853396', c(0:5))
samples_DRG_single <- paste0('SRR85339', c(66:86))

samples_DRG_paired
samples_DRG_single

# trimming - single-----
dir <- '/data/NCATS_ifx/data/dbGaP/DRGs/SRA/FASTQs'
samples <- unique(samples_DRG_single)
for (sample in samples){
  cat(paste0("java -jar $TRIMMOJAR SE -phred33 ", dir ,'/', sample, '.sra.fastq ', sample, '.trim.fastq.gz ',
             'ILLUMINACLIP:/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 LEADING:10 TRAILING:5 MAXINFO:50:0.97 MINLEN:36'), sep="\n")
}

# trimming - paired -----
dir <- '/data/NCATS_ifx/data/dbGaP/DRGs/SRA/FASTQs'
samples <- unique(samples_DRG_paired)
for (sample in samples){
  cat(paste0("java -jar $TRIMMOJAR PE -phred33 ", dir ,'/', sample, '.sra_1.fastq ', dir ,'/', sample, '.sra_2.fastq -baseout ', sample, '.trim.fastq.gz ',
             'ILLUMINACLIP:/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:5 MAXINFO:50:0.97 MINLEN:36'), sep="\n")
}

# swarm -f trim_swarm_paired.sh -g 6 -t 8 --module trimmomatic --time=24:00:00
# swarm -f trim_swarm_single.sh -g 6 -t 8 --module trimmomatic --time=24:00:00

# align, single----
indir <- '/data/NCATS_ifx/data/dbGaP/DRGs/SRA/Trimmed_FASTQs'
outdir <- '/data/NCATS_ifx/data/dbGaP/DRGs/SRA/Aligned'
samples <- unique(samples_DRG_single)
for (sample in samples){
  cat(paste0('cd ', outdir, ' && STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /fdb/STAR_current/GENCODE/Gencode_human/release_27/genes-75 --sjdbOverhang 75 --outSAMunmapped Within --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesIn ', indir, '/' , sample , '.trim.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix Sample_' , sample , '_hg38'))
  cat("\n\n")
}

#
# swarm -f star_single.sh -t 12 -g 40 --module STAR --time=24:00:00 # running on biowulf.

# align, paired----
samples <- unique(samples_DRG_paired)
for (sample in samples){
  cat(paste0('cd ', outdir, ' && STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /fdb/STAR_current/GENCODE/Gencode_human/release_27/genes-75 --sjdbOverhang 75 --outSAMunmapped Within --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesIn ', indir, '/' , sample , '.trim_1P.fastq.gz ', indir, '/' , sample , '.trim_2P.fastq.gz',' --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix Sample_' , sample , '_hg38'))
  cat("\n\n")
}

# --readFilesIn ', dir, '/' , sample , '_1P.fastq.gz ', dir, '/', sample , '_2P.fastq.gz

# swarm -f star_paired.sh -t 12 -g 40 --module STAR --time=24:00:00

# htseq-count----
indir <- '/data/NCATS_ifx/data/dbGaP/DRGs/SRA/Aligned/'
outdir <- '/data/NCATS_ifx/data/dbGaP/DRGs/SRA/htseq'

samples <- c(samples_DRG_paired,  samples_DRG_single)
samples <- samples[samples %nin% c('SRR8533978','SRR8533976','SRR8533972')]

for (sample in samples){
  cat(paste0('htseq-count -f bam -r pos -s no -t exon -m union ', indir , 'Sample_' , sample , '_hg38Aligned.sortedByCoord.out.bam /fdb/GENCODE/Gencode_human/release_27/gencode.v27.annotation.gtf > ', outdir, '/', sample , '_htseq_counts.txt' , "\n\n"))
}
#swarm -f htseq.sh -t 4 -g 4 --module htseq --time=24:00:00 #


#######################

# convert ENSG to gene symbol
library(biomaRt)
library(stringr)
SRR_ex <- fread('~/Downloads/DRGs/Countfiles_ENSG/SRR8533960_htseq_counts.txt', header=F)

SRR_ex <- SRR_ex[1:(nrow(SRR_ex)-5),]
SRR_ex
dt <- SRR_ex
names(dt) <- c("ENSG.full", "Counts")

ensg.genes <- data.table("ENSG.full" = dt$ENSG.full)

ensg.genes[,ENSG.short := tstrsplit(ENSG.full, "\\.")[1]]

genes <- ensg.genes$ENSG.short
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),
                values=genes,mart= mart)
G_list <- as.data.table(G_list)
fwrite(G_list, '/Volumes/ncatssctl/NGS_related/Genomes/Biomart_hsapiens_ensembl_gene_symbols.csv')
fwrite(G_list, '~/Downloads/DRGs/Biomart_hsapiens_ensembl_gene_symbols.csv')
ensg.genes <- merge(ensg.genes, G_list, all=T, by.x="ENSG.short", by.y="ensembl_gene_id")
ensg.genes <- na.omit(ensg.genes)
dt <- subset(dt, dt$ENSG.full %in% ensg.genes$ENSG.full)
dt <- merge(dt, ensg.genes, by="ENSG.full")
dt <- dt[,c("external_gene_name", "Counts")]
dt <- dt[order(external_gene_name, -Counts)]
dt <- dt[!duplicated(external_gene_name),]

sample_id <- 'SRR8533960'
fwrite(dt, paste0("~/Downloads/DRGs/Countfiles_gene_symbol/", sample_id, "_gene_symbol_counts.txt"), col.names = F, row.names = F, quote=F, sep="\t")
dt <- fread('~/Downloads/DRGs/Countfiles_gene_symbol/SRR8533960_gene_symbol_counts.txt')

setwd('~/Downloads/DRGs/Countfiles_ENSG/')
files <- Sys.glob('*_htseq_counts.txt')
files.table <- data.table(file=files)
files.table[,Run:= tstrsplit(file,'_')[1]]
files.table

for (file in files.table$file){
  dt <- fread(file)
  sample_id <- str_split_fixed(file, "_htseq_counts.txt",2)[1]
  names(dt) <- c("ENSG.full", "Counts")
  dt <- dt[ENSG.full %in% ensg.genes$ENSG.full,]
  dt[,external_gene_name := ensg.genes$external_gene_name]
  dt <- subset(dt, select=c('external_gene_name', 'Counts'))
  dt <- dt[order(external_gene_name, -Counts)]
  dt <- dt[!duplicated(external_gene_name),]
  fwrite(dt, paste0('~/Downloads/DRGs/Countfiles_gene_symbol/',sample_id,"_gene_symbol_counts.txt"), col.names = F, row.names = F, quote=F, sep="\t")
}

TaoNoci_ex <- fread('~/Downloads/DRGs/ISB012/H9noc_D0_1_gene_symbol_counts.txt', header=F)

intersection <- intersect(TaoNoci_ex$V1, dt$V1)
intersection

setwd('../Countfiles_gene_symbol/')
files.table <- data.table(file=Sys.glob('*txt'))
files.table[,Run:= tstrsplit(file,'_')[1]]

TaoNoci_ex <- subset(TaoNoci_ex, TaoNoci_ex$V1 %in% intersection)
dt <- fread('~/Downloads/DRGs/Countfiles_gene_symbol/SRR8533960_gene_symbol_counts.txt')
names(dt) <- c('external_gene_name', 'Counts')
dt <- subset(dt, select=c('external_gene_name', 'Counts'))
dt <- dt[order(external_gene_name, -Counts)]
dt <- dt[!duplicated(external_gene_name),]
dt <- subset(dt, dt$external_gene_name %in% intersection)
all(dt$V1 == TaoNoci_ex$V1)

for (file in files.table$file){
  dt <- fread(file, header=F)
  sample_id <- str_split_fixed(file, "_",Inf)[1]
  names(dt) <- c('external_gene_name', 'Counts')
  dt <- dt[order(external_gene_name, -Counts)]
  dt <- dt[!duplicated(external_gene_name),]
  dt <- subset(dt, dt$external_gene_name %in% intersection)
  fwrite(dt, paste0('~/Downloads/DRGs/Countfiles_gene_symbol/intersection/', sample_id,"_intersection_gene_symbol_counts.txt"), col.names = F, row.names = F, quote=F, sep="\t")
}

setwd('../ISB012/')
files.table <- data.table(file=Sys.glob('*txt'))
files.table[,Sample:= tstrsplit(file,'_gene_symbol_counts.txt')[1]]
files.table

for (file in files.table$file){
  dt <- fread(file, header=F)
  sample_id <- str_split_fixed(file, "_gene_symbol_counts.txt",Inf)[1]
  names(dt) <- c('external_gene_name', 'Counts')
  dt <- dt[order(external_gene_name, -Counts)]
  dt <- dt[!duplicated(external_gene_name),]
  dt <- subset(dt, dt$external_gene_name %in% intersection)
  fwrite(dt, paste0('~/Downloads/DRGs/Countfiles_gene_symbol/intersection/', sample_id,"_intersection_gene_symbol_counts.txt"), col.names = F, row.names = F, quote=F, sep="\t")
}

isb012.ex <- fread('~/Downloads/DRGs/Countfiles_gene_symbol/intersection/H9noc_D0_1_intersection_gene_symbol_counts.txt', header=F)
isb012.ex

drg.ex <- fread('~/Downloads/DRGs/Countfiles_gene_symbol/intersection/SRR8533960_intersection_gene_symbol_counts.txt', header=F)
drg.ex
all(drg.ex$V1 == isb012.ex$V1)
all(isb012.ex$V1 == drg.ex$V1)
all(isb012.ex$V1 == intersection)
all(drg.ex$V1 == intersection)

length(drg.ex$V1)
length(isb012.ex$V1)

# merge into one table because deseq is complaining.
setwd('../Countfiles_gene_symbol/intersection/')

counts <- data.table(GeneId = intersection)

for (file in Sys.glob('*txt')){
  counts.tmp <- fread(file, header=F)
  sample_id <- str_split_fixed(file, '_intersection_gene_symbol_counts.txt', Inf)[1]
  names(counts.tmp)<- c('GeneId', sample_id)
  counts <- merge(counts, counts.tmp, by='GeneId', all.x=T, all.y=F)
}

cat(names(counts), sep='\n')
sampleTable <- as.data.table(read_xlsx('~/Downloads/DRGs/DRGs_merged_sampletable.xlsx'))
sampleTable
colData <- data.frame(condition=sampleTable$condition, row.names = sampleTable$sample_id)
counts <- as.data.frame(counts)
row.names(counts) <- counts$GeneId
counts <- counts[-1]
counts[1:10,1:10]
colData
names(counts)
counts <- subset(counts, select=c(1:3, 13:18, 4:12, 19:42))
counts[1:10,1:10]
counts <- na.omit(counts)

dds <- DESeqDataSetFromMatrix(countData = counts,colData = colData, design = ~ condition)
dds
