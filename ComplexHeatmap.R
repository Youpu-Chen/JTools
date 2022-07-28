rm(list = ls())
options(stringsAsFactors = F)

library(DESeq2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

# Load data
rawdata <- read.csv('pandas_all_samples.csv', header = TRUE, row.names = 1)
rawdata[is.na(rawdata)] <- 0

# Get the sample name
names(rawdata)
table(table(names(rawdata)))  # count sample number
ncol(rawdata)

# Create the design matrix
designmatrix = data.frame(names(rawdata), factor(rep(c('RT', 'ST', 'R', 'T'), each = 3)))
names(designmatrix) <- c('sample', 'group')
nrow(designmatrix)

# Check whether all the samples of rawdata are in the designmatrix
all(colnames(rawdata) == designmatrix$sample)
all(colnames(rawdata) %in% designmatrix$sample)


# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=as.matrix(rawdata), 
                              colData=designmatrix, 
                              design=~group)

# Pre-filtering
dim(dds[rowSums(counts(dds)) > 5, ])
dds <- dds[rowSums(counts(dds)) > 5, ]  # filtering the gene whose total read counts of all samples are smaller than 5

# Run DESeq2
dds <- DESeq(dds)

# Get results info
res <- results(dds)
head(results(dds))

# Summary of differential gene expression
summary(res)

# Sort summary list by p-value
resOrdered <- res[order(res$padj),]

# Save significantly expressed gene
sig <- resOrdered[!is.na(resOrdered$padj) &
    resOrdered$padj < 0.05 &
    abs(resOrdered$log2FoldChange)>=2, ]
# write.csv(sigs, file = "deseq_results.csv")

# Create a log2FoldChange dataset
l2_val <- as.matrix(sig$log2FoldChange)
colnames(l2_val) <- "logFC"

# Plot heatmap
h <- Heatmap(l2_val, # row_labels = df.top$symbol[rows_keep], 
            cluster_rows = F, name="logFC", # top_annotation = ha, 
            # col = col_logFC,
            # cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            # grid.text(round(l2_val[i, j],2), x, y)
            # }
            )
h