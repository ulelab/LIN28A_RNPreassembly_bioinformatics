library( "DESeq2" )
library("ggplot2")
library("apeglm")
library("tximport")
library("stringr")

# Import quantfiles
quantfiles <- c(Sys.glob('../../data/MihaDeseq/salmon_quantfiles/KO*FCL*'), Sys.glob('../../data/MihaDeseq/salmon_quantfiles/S200*FCL*'))
names <- str_split_i(quantfiles, pattern = '/', i = -1) %>% str_split_i(pattern = '.quant', i = 1)
# Set names for quantfiles
names(quantfiles) <- names
print(quantfiles)
# print(names)


#Import transcript to gene mappings for tximport
tx_gene <- read.csv('transcript_to_gene_mapping_mm10_V22.csv', header = TRUE, sep = ",")
head(tx_gene)

count_data <- tximport(files = quantfiles,
         type = "salmon",
         tx2gene = tx_gene,
         ignoreTxVersion = TRUE
         )

countData <- as.data.frame(count_data[["counts"]])
countData <- round(countData, digits = 0)
# Make a metadata table with matching sample names
genotypes <- str_split_i(colnames(countData), pattern = '_', i = 1)
metaData <- data.frame(row.names = colnames(countData), genotype = genotypes)
metaData$genotype <- factor(metaData$genotype)
print(metaData)


print("Creating deseq2 object")
dds <- DESeqDataSetFromMatrix(countData=countData,
                              colData=metaData,
                              design=~genotype, tidy = FALSE)
dds <- estimateSizeFactors(dds)
print(dds)

normalized_counts <- counts(dds, normalized=TRUE)

# For each condition get rows where at lest 2 replicates have a norm count > 5
# keep and idx gave the same result

idx <- c()
for (g in unique(metaData$genotype))
{
  print(g)
  # Get the names of rows with the same genotype
  cols <- rownames(metaData)[which(metaData$genotype == g)]
  # Get dataframe with relevant columns
  dfTemp <- normalized_counts[, cols]
  # Apply filter - get rows where value is greater than 5 in at least 2 columns
  # keepRows <- rownames(dfTemp)[which(rowSums( dfTemp >= 5 ) >= 2)]
  id <- which(rowSums( dfTemp >= 5 ) >= 2)
  # keep <- c(keep, keepRows)
  idx <- c(idx, id)
  # print(cols)
}
keep <- unique(idx)
# Filter dds
dds <- dds[keep,]
print(dds)

## Because we have quantseq data, we don't need the length adjustment with Tximport (would be important for normal RNAseq)
# dds <- DESeqDataSetFromTximport(count_data, metaData, ~genotype)

dds <- DESeq(dds)

# Save normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, file='normalized_counts.csv', row.names=TRUE)

# Save results
res_KO_WT <- lfcShrink(dds, coef="genotype_S200WT_vs_KO", type="apeglm")
# res_KO_WT <- results(dds, contrast=c("genotype","KO",  "S200WT"))
res_KO_WT <- res_KO_WT[order(res_KO_WT$padj),]
write.csv(res_KO_WT, file='res_S200WT_vs_KO_deseq_lfcShrink.csv', row.names=TRUE)

res_KO_S200A <- lfcShrink(dds, coef="genotype_S200A_vs_KO", type="apeglm")
# res_KO_S200A <- results(dds, contrast=c("genotype","KO",  "S200A"))
res_KO_S200A <- res_KO_S200A[order(res_KO_S200A$padj),]
write.csv(res_KO_S200A, file='res_S200A_vs_KO_deseq_lfcShrink.csv', row.names=TRUE)


head(res_KO_WT) #let's look at the results table

summary(res_KO_WT) #summary of results


# Make a basic volcano plot

# KO vs WT
pdf('volcano_KO_vs_S200WT.pdf')
with(res_KO_WT, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot KO vs S200WT"))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res_KO_WT, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res_KO_WT, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

# KO vs S200A
pdf('volcano_KO_vs_S200A.pdf')
with(res_KO_S200A, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot KO vs S200A"))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res_KO_S200A, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res_KO_S200A, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

# PCA plot

#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation
vsdata <- vst(dds, blind=FALSE)
# pdf('PCA_afterVST.pdf', width=5, height=5)
p <- plotPCA(vsdata, intgroup="genotype")
p <- p + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black", position = position_dodge(width = 1))+ coord_fixed()
print(p)
ggsave('PCA_afterVST.pdf', height=3, width=4, scale=2, dpi=300)
# dev.off()

# rlog transform
rld <- rlog(dds)
# pdf('PCA_afterRlogTransform.pdf', width=5, height=5)
p <- plotPCA(rld, intgroup="genotype")
p <- p + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black", position = position_dodge(width = 1)) + coord_fixed()
print(p)
ggsave('PCA_afterRlogTransform.pdf', height=3, width=4, scale=2, dpi=300)
# dev.off()
