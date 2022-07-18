# This script is used to conduct the differential expression analysis
rm(list = ls())
options(stringsAsFactors = F)

#-----------------------------------
# Preparation
#-----------------------------------
# Load packages
# DESeq links: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
# BiocManager::install("DESeq2")
library(DESeq2)

# Load data
# Note: if your sample name start with number, R will add X in front of it,
# because using number-started string as variable name is not allowed in R
rawdata <- read.csv('pandas_all_samples.csv', header = TRUE, row.names = 1)
# if your raw count have NA in the dataframe, please run the code below,
# if it is not, you do not have to run it.
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


#-----------------------------------
# DESeq2 analysis
#-----------------------------------
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
res <- res[order(res$padj),]
head(res)


#-----------------------------------
# Visualization
#-----------------------------------
library(ggplot2)
library(scales) # needed for oob parameter
library(viridis)

# MA plot
deseq2ResDF <- as.data.frame(res)
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA)
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) + 
	geom_point(size=1) + 
	scale_y_continuous(limits=c(-3, 3), oob=squish) + 
	scale_x_log10() + 
	geom_hline(yintercept = 0, colour="tomato1", size=2) + 
	labs(x="mean of normalized counts", y="log fold change") + 
	scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") +
	theme_bw()

# Viewing normalized counts for a single geneID
# Note: code below could be used to plot specified gene expression in all the groups
otop2Counts <- plotCounts(dds, gene="osa-miR1861b", intgroup=c("sample", "group"), returnData=TRUE)

colourPallette <- c("#7145cd","#bbcfc4","#90de4a","#cd46c1","#77dd8e","#592b79","#d7c847","#6378c9","#619a3c","#d44473","#63cfb6","#dd5d36")

ggplot(otop2Counts, aes(x=group, y=count, colour=sample, group=sample)) + 
	geom_point() + 
    # Note: when you have multiple tissues from one sample, you can uncomment the line below
	# geom_line() +  
	theme_bw() + 
	theme(axis.text.x=element_text(angle=15, hjust=1)) + 
	scale_colour_manual(values=colourPallette) +
	guides(colour=guide_legend(ncol=3)) + 
	ggtitle("OTOP2")

# Visualizing expression data with a heatmap
# Transform count data using the variance stablilizing transform
deseq2VST <- varianceStabilizingTransformation(dds)

# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)

# Keep only the significantly differentiated genes where the fold-change was at least 3
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05 & abs(deseq2ResDF$log2FoldChange) > 3,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

# Convert the VST counts to long format for ggplot2
library(reshape2)

# First compare wide vs long version
deseq2VST_wide <- deseq2VST
deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))

head(deseq2VST_wide)
head(deseq2VST_long)

# Now overwrite our original data frame with the long format
deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))

# Make a heatmap
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + 
	geom_raster() + 
	scale_fill_viridis(trans="sqrt") + 		
	theme(axis.text.x=element_text(angle=65, hjust=1), 
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
heatmap

# Ref
# [1] https://lashlock.github.io/compbio/R_presentation.html
# [2] https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/