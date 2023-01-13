# install if not already installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

# load the library
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene


# obtain gene locations from the TxDb object as a data frame
genes <- genes(TxDb, columns=c("gene_id"), single.strand.genes.only=FALSE)
genes <- as.data.frame(genes)
head(genes)

chromosomes <- paste0("chr", c(1:22,"X","Y"))
genes <- genes[genes$seqnames %in% chromosomes,]
genes <- unique(genes[,c("group_name", "seqnames", "start", "end", "width", "strand")])
colnames(genes) <- c("gene_id", "seqnames", "start", "end", "width", "strand")
head(genes)
table(genes$gene_id)


immuneGenes <- read.delim("http://genomedata.org/gen-viz-workshop/intro_to_ggplot2/immuneGenes.tsv")

# merge the immune gene annotations onto the gene object
genes <- merge(genes, immuneGenes, by.x=c("gene_id"), by.y=c("entrez"), all.x=TRUE)
genesOnChr6 <- genes[genes$seqnames %in% "chr6",]
head(genesOnChr6)


ggplot(data=genes, aes(x=start + .5*width,)) + 
    geom_density() + 
    xlab("Chromosomal position of genes (gene start + 0.5 * width)")

