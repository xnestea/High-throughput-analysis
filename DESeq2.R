library(DESeq2)
setwd("Desktop/Mouse_organ_development/")
counts <- as.matrix(read.csv(file = "genes.raw.htseq2.tsv",sep = "\t", row.names = 1))
metadata <- read.csv("E-MTAB-2328.sdrf.txt", sep = "\t")
metadata <- subset(metadata, Characteristics.developmental.stage. %in% c("embryonic day 15.5", "postnatal day 29"))
metadata <- subset(metadata, Characteristics.organism.part. %in% c("liver"))
list_of_names <- metadata$Scan.Name
list_of_names <- substr(list_of_names, start = 1, stop = 6)
counts <- subset(counts, colnames(counts) %in% c(list_of_names))
#deseq <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design= ~Characteristics.developmental.stage.)

