# ------------------------------------------------------------------------------
# Description: Transform single sample read counts to RPKM and filtering analysis
# for GO enrichment analysis
# Date: 2022-12-18
# ------------------------------------------------------------------------------
rm(list=ls())
options(stringsAsFactors = F)

# load data
counts <- read.delim('rna.mat', comment.char = '#', header = TRUE, sep = '\t')
names(counts)[7] <- "S1"
head(counts)

# RPKM
total.reads <- sum(counts[, 7])
rpkm.line <- function(line){
  return( as.numeric(line[7])/( as.numeric(line[6])/1000 * total.reads/1000000) )
}

counts$rpkm <- apply(counts, MARGIN = 1, FUN = rpkm.line)
head(counts)

# retrieve the top 500 level of genes
counts <- counts[order(counts$rpkm, decreasing = TRUE), ]
head(counts)
counts$Geneid[1:500]
if(F){
  test.vector <- as.vector(counts$Chr[1:500])
  selectOne <- function(value) {
    strsplit(value, split = ';')[[1]][1]
  }
  test.df <- data.frame(table(unlist(lapply(test.vector, FUN = selectOne))))
  test.df <- test.df[order(test.df$Freq, decreasing = TRUE), ]
  head(test.df)
  # Var1 Freq
  # 4    chrIV   74
  # 12  chrXII   58
  # 8   chrVII   52
  
  strsplit(test.vector[450], split = ';')[[1]][1]
}


# enrichment analysis
# BiocManager::install("org.Sc.sgd.db")
library(org.Sc.sgd.db)
library(clusterProfiler)
library(ggplot2)

# GO
ego <- enrichGO(gene = counts$Geneid[1:500],
                OrgDb = org.Sc.sgd.db,
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 10,
                pvalueCutoff = 0.05,
                keyType='ENSEMBL')

p1 <- barplot(ego, showCategory = 10)
p1
p2 <- dotplot(ego, showCategory = 10)
p2
# ggsave("GO-enrich-dotplot.pdf", p2, height = 6, width = 10)


# KEGG
ekegg <- enrichKEGG(gene = counts$Geneid[1:500],
                    organism = 'sce',
                    pAdjustMethod = "BH",
                    minGSSize = 10,
                    pvalueCutoff = 0.05,
                    keyType='kegg')

p3 <- emapplot(ekegg, showCategory = 20)
p3
# ggsave("KEGG-enrich-network.pdf", p3, height = 6, width = 6)



### top3 GO functions
library(tidyverse)
library(DOSE)
library(GOSemSim)
if(F){
  # required functions
  get_GO_data <- function(OrgDb, ont, keytype) {
    GO_Env <- get_GO_Env()
    use_cached <- FALSE
    
    if (exists("organism", envir=GO_Env, inherits=FALSE) &&
        exists("keytype", envir=GO_Env, inherits=FALSE)) {
      
      org <- get("organism", envir=GO_Env)
      kt <- get("keytype", envir=GO_Env)
      
      if (org == DOSE:::get_organism(OrgDb) &&
          keytype == kt &&
          exists("goAnno", envir=GO_Env, inherits=FALSE)) {
        ## https://github.com/GuangchuangYu/clusterProfiler/issues/182
        ## && exists("GO2TERM", envir=GO_Env, inherits=FALSE)){
        
        use_cached <- TRUE
      }
    }
    
    if (use_cached) {
      goAnno <- get("goAnno", envir=GO_Env)
    } else {
      OrgDb <- GOSemSim:::load_OrgDb(OrgDb)
      kt <- keytypes(OrgDb)
      if (! keytype %in% kt) {
        stop("keytype is not supported...")
      }
      
      kk <- keys(OrgDb, keytype=keytype)
      goAnno <- suppressMessages(
        AnnotationDbi::select(OrgDb, keys=kk, keytype=keytype,
                              columns=c("GOALL", "ONTOLOGYALL")))
      
      goAnno <- unique(goAnno[!is.na(goAnno$GOALL), ])
      
      assign("goAnno", goAnno, envir=GO_Env)
      assign("keytype", keytype, envir=GO_Env)
      assign("organism", DOSE:::get_organism(OrgDb), envir=GO_Env)
    }
    
    if (ont == "ALL") {
      GO2GENE <- unique(goAnno[, c(2,1)])
    } else {
      GO2GENE <- unique(goAnno[goAnno$ONTOLOGYALL == ont, c(2,1)])
    }
    
    GO_DATA <- DOSE:::build_Anno(GO2GENE, get_GO2TERM_table())
    
    goOnt.df <- goAnno[, c("GOALL", "ONTOLOGYALL")] %>% unique
    goOnt <- goOnt.df[,2]
    names(goOnt) <- goOnt.df[,1]
    assign("GO2ONT", goOnt, envir=GO_DATA)
    return(GO_DATA)
  }
  
  get_GO_Env <- function () {
    if (!exists(".GO_clusterProfiler_Env", envir = .GlobalEnv)) {
      pos <- 1
      envir <- as.environment(pos)
      assign(".GO_clusterProfiler_Env", new.env(), envir=envir)
    }
    get(".GO_clusterProfiler_Env", envir = .GlobalEnv)
  }
  
  get_GO2TERM_table <- function() {
    GOTERM.df <- get_GOTERM()
    GOTERM.df[, c("go_id", "Term")] %>% unique
  }
  
  get_GOTERM <- function() {
    pos <- 1
    envir <- as.environment(pos)
    if (!exists(".GOTERM_Env", envir=envir)) {
      assign(".GOTERM_Env", new.env(), envir)
    }
    GOTERM_Env <- get(".GOTERM_Env", envir = envir)
    if (exists("GOTERM.df", envir = GOTERM_Env)) {
      GOTERM.df <- get("GOTERM.df", envir=GOTERM_Env)
    } else {
      GOTERM.df <- toTable(GOTERM)
      assign("GOTERM.df", GOTERM.df, envir = GOTERM_Env)
    }
    return(GOTERM.df)
  }
  
}

GO.BP.dat <- get_GO_data("org.Sc.sgd.db", "BP", "ENSEMBL")

### peptide metabolic process
BP.peptide_metabolic_process <- names(GO.BP.dat$PATHID2NAME[grep("peptide metabolic process", GO.BP.dat$PATHID2NAME)])
ENSEMBL.BP.peptide_metabolic_process <- unique(unlist(GO.BP.dat$PATHID2EXTID[ BP.peptide_metabolic_process ]))

# counts$Geneid[1:500] %in% ENSEMBL.BP.peptide_metabolic_process
geneset1 <- counts[1:500, c(1,2)][ counts$Geneid[1:500] %in% ENSEMBL.BP.peptide_metabolic_process, ]
# length(geneset1)

### peptide biosynthetic process
BP.peptide_biosynthetic_process <- names(GO.BP.dat$PATHID2NAME[grep("peptide biosynthetic process", GO.BP.dat$PATHID2NAME)])
ENSEMBL.BP.peptide_biosynthetic_process <- unique(unlist(GO.BP.dat$PATHID2EXTID[ BP.peptide_biosynthetic_process ]))

# counts$Geneid[1:500] %in% ENSEMBL.BP.peptide_biosynthetic_process
geneset2 <- counts[1:500, c(1,2)][ counts$Geneid[1:500] %in% ENSEMBL.BP.peptide_biosynthetic_process, ]
# length(geneset2)


### translation
BP.translation <- names(GO.BP.dat$PATHID2NAME[grep("^translation", GO.BP.dat$PATHID2NAME)])
ENSEMBL.BP.translation <- unique(unlist(GO.BP.dat$PATHID2EXTID[ BP.translation ]))

geneset3 <- counts[1:500, c(1,2)][ counts$Geneid[1:500] %in% ENSEMBL.BP.translation, ]
head(geneset3)

# acid
BP.acid <- names(GO.BP.dat$PATHID2NAME[grep("carboxylic acid biosynthetic process", GO.BP.dat$PATHID2NAME)])
ENSEMBL.BP.acid <- unique(unlist(GO.BP.dat$PATHID2EXTID[ BP.acid ]))

geneset.acid <- counts[1:500, c(1,2)][ counts$Geneid[1:500] %in% ENSEMBL.BP.acid, ]
head(geneset.acid)
geneset.acid
write.csv(geneset.acid, 'acid-related-gene.csv', quote = F, row.names = F)

#
which.max(table(unlist(lapply(geneset1$Chr, FUN = selectOne))))
# chrIV
# chrXII
# chrVII

which.max(table(unlist(lapply(geneset2$Chr, FUN = selectOne))))
# chrIV 
# chrXII
# chrVII

which.max(table(unlist(lapply(geneset3$Chr, FUN = selectOne))))
# chrIV 
# chrXII
# chrVII


### core gene set
# 
geneset1.chrIV <- geneset1[(geneset1$Chr == "chrIV" | geneset1$Chr == "chrIV;chrIV"), 1]
geneset2.chrIV <- geneset2[(geneset2$Chr == "chrIV" | geneset2$Chr == "chrIV;chrIV"), 1]
geneset3.chrIV <- geneset3[(geneset3$Chr == "chrIV" | geneset3$Chr == "chrIV;chrIV"), 1]
chrIV.name <- unique(c(geneset1.chrIV, geneset2.chrIV, geneset3.chrIV))

# 
geneset1.chrXII <- geneset1[(geneset1$Chr == "chrXII" | 
                               geneset1$Chr == "chrXII;chrXII" | 
                               geneset1$Chr == "chrXII;chrXII;chrXII"), 1]
geneset2.chrXII <- geneset2[(geneset2$Chr == "chrXII" | 
                               geneset2$Chr == "chrXII;chrXII" | 
                               geneset2$Chr == "chrXII;chrXII;chrXII"), 1]
geneset3.chrXII <- geneset3[(geneset3$Chr == "chrXII" | 
                               geneset3$Chr == "chrXII;chrXII" | 
                               geneset3$Chr == "chrXII;chrXII;chrXII"), 1]

chrXII.name <- unique(c(geneset1.chrXII, geneset2.chrXII, geneset3.chrXII))

# 
geneset1.chrVII <- geneset1[(geneset1$Chr == "chrVII" | 
                               geneset1$Chr == "chrVII;chrVII" | 
                               geneset1$Chr == "chrVII;chrVII;chrVII"), 1]
geneset2.chrVII <- geneset2[(geneset2$Chr == "chrVII" | 
                               geneset2$Chr == "chrVII;chrVII" | 
                               geneset2$Chr == "chrVII;chrVII;chrVII"), 1]
geneset3.chrVII <- geneset3[(geneset3$Chr == "chrVII" | 
                               geneset3$Chr == "chrVII;chrVII" | 
                               geneset3$Chr == "chrVII;chrVII;chrVII"), 1]
chrVII.name <- unique(c(geneset1.chrVII, geneset2.chrVII, geneset3.chrVII))


### test
which.min(aggregate(counts$rpkm~counts$Chr, FUN =sum))
