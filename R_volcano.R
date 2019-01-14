if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggplot2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("gridExtra")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("RColorBrewer")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("reshape2")


library(readxl)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(reshape2)

setwd("c:/Users/xNesTea/Documents/project_HT/Scripts/High-throughput-analysis/")
DE_table <- as.data.frame(read_excel("./DESEQ2_volcano.xlsx"))
rownames(DE_table) = DE_table$Gene
DE_table$threshold = as.factor(DE_table$padj < 0.05)
##Construct the plot object
g <- ggplot(data=DE_table, aes(x=DE_table$log2FoldChange, y=-log10(DE_table$padj)) )+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,color='gray16'),
        axis.title.x = element_text(color='gray45',size=15,family="Times",margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(color='gray45',size=15,family="Times",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(color='gray45',size=15,family="Times"),
        legend.title=element_text(color='gray45',size=15,family="Times"),
        legend.text=element_text(color='gray45',size=15,family="Times"),
        legend.position="none") +
  geom_point(alpha=0.7, size=1.75,shape=16, aes(colour = DE_table$threshold)) +
  scale_color_manual(values = c("grey", "coral")) +
  #geom_text(aes(label = ifelse(gene_list3$threshold == TRUE, rownames(gene_list3),"")),size=3,check_overlap = TRUE,hjust = 0, nudge_x = 0.05) +
  xlim(c(-3, 3)) + ylim(c(0, 7.5)) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted P-value")  +
  ggtitle("DE_volcano")

ggsave("./DE_plot", plot = g)