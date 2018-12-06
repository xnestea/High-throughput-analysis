### install necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")

### load library
library(DESeq2)

### set working directory to the folder containg data files (matadata and so on)
setwd("project_HT/Data/")

### read table with number of reads for genes in form of matrix
counts <- as.matrix(read.csv(file = "raw/genes.raw.htseq2.tsv",sep = "\t", row.names = 1))

### read in metadata file and restrict it to data we are interested in
### also drop technical duplicates since there are not present in count matrix
metadata <- read.csv("E-MTAB-2328.sdrf.txt", sep = "\t")
metadata <- subset(metadata, Characteristics.developmental.stage. %in% c("embryonic day 15.5", "postnatal day 29"))
metadata <- subset(metadata, Characteristics.organism.part. %in% c("liver"))
metadata <- metadata[!duplicated(metadata$Assay.Name),]
### get sample names we are interested in from restricted metadata
### and restrict our read count matrix t only have samples we are interested in
list_of_names <- metadata$Scan.Name
list_of_names <- substr(list_of_names, start = 1, stop = 6)
counts <- counts[, colnames = colnames(counts) %in% c(list_of_names)]

#to be continued

### perform DESeq analysis
deseq <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design= ~Characteristics.developmental.stage.)
analysis <- DESeq(deseq)
### visualizing results to be clearly visible.
res = results(analysis, 
              contrast=c("Characteristics.developmental.stage.","embryonic day 15.5","postnatal day 29"), 
              lfcThreshold = 0,altHypothesis="greaterAbs",
              alpha = 0.01)

write.table(res ,file = "DESEQ2_out.txt", sep = "\t",row.names = T,col.names = NA)
### PIANO PART


