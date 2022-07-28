# Load packages
library(pheatmap)
library(RColorBrewer)


if(F){
    # Load raw Count Matrix
    rawdata <- read.csv('pandas_all_samples.csv', header = TRUE, row.names = 1)
    rawdata[is.na(rawdata)] <- 0
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
    res <- results(dds) # contrast=c('condition') should be specified

    # Sort summary list by p-value
    resOrdered <- res[order(res$padj),]
    head(res)

    # Save significantly expressed gene
    sig <- resOrdered[!is.na(resOrdered$padj) &
        resOrdered$padj < 0.05 &
        abs(resOrdered$log2FoldChange)>=2, ]
    dim(sig)
    head(sig)
}


# Using significant expressed genes as the input of heatmap
sigene <- rownames(sig)
sigmatrix <- rawdata[sigene,]

# Using log2
pheatmap(sigmatrix)
# pheatmap(log2(sigmatrix))
