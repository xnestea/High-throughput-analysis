library(biomaRt)
### table for mapping ensembl to gene symbol for our species - Mus Musculus
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

map_gene <- getBM(attributes=c('ensembl_gene_id',
                          'external_gene_name'),
             mart = ensembl)

### formating the table for later use
colnames(map_gene)[colnames(map_gene)=="external_gene_name"] <- "Gene"
map_gene$Gene <- toupper(map_gene$Gene)
rownames(map_gene) <- map_gene[,1]
map_gene[,1] <- NULL

library(piano)
### read in the file from DESeq2
DE <- read.delim(file="DESEQ2_out.txt", row.names = 1, stringsAsFactors = F)

### make a new row in DE table for respective Gene symbol 
DE[,"Gene"] <- map_gene[row.names(DE), "Gene"]

### remove rows for which no gene symbols is available
DE <- DE[!is.na(DE$Gene),]

### remove rows with same gene symbols
DE <- DE[!duplicated(DE$Gene),]

###Gene symbols are now row names, replacing Ensembl entries
row.names(DE) <- DE$Gene

### we only need logfoldchange and pvalue as piano input
DE <- DE[ ,c('log2FoldChange','pvalue')] 

### make both logdfoldchange and p values as separate matrices
### with Genesymbols and rownames
pval <- as.matrix(DE[ ,2])
lfc <- as.matrix(DE[ ,1])
row.names(pval) <- row.names(DE)
row.names(lfc) <- row.names(DE)

### substitute NA with reasonable values for piano
pval[is.na (pval)] <- 1
lfc[is.na (lfc)] <- 0

### load gene set enrichement grouping, downloaded from go2msig.org
gset=loadGSC("../../Data/Mus_musculus_GSEA_GO_sets_all_symbols_April_2015.gmt") 

### this is GSE analysis, adjMethod flag means we want to get fdr
gsaRes <- runGSA(pval,lfc,gsc=gset, nPerm = 1000, adjMethod = "fdr")

### Just a simple write function for GSE output
GSAsummaryTable(gsaRes, save=TRUE, file="piano.txt")