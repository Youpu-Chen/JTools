rm(list=ls())
options(stringsAsFactors = F)


#install if necessary
# source("http://bioconductor.org/biocLite.R")
# biocLite("Rsamtools")

#load library
library(Rsamtools)

#read in entire BAM file
bam <- scanBam("wgEncodeRikenCageHchCellPapAlnRep1.bam")

#names of the BAM fields
names(bam[[1]])
# [1] "qname"  "flag"   "rname"  "strand" "pos"    "qwidth" "mapq"   "cigar"
# [9] "mrnm"   "mpos"   "isize"  "seq"    "qual"

#distribution of BAM flags
table(bam[[1]]$flag)

#      0       4      16 
#1472261  775200 1652949

#function for collapsing the list of lists into a single list
#as per the Rsamtools vignette
.unlist <- function (x){
   ## do.call(c, ...) coerces factor to integer, which is undesired
   x1 <- x[[1L]]
   if (is.factor(x1)){
      structure(unlist(x), class = "factor", levels = levels(x1))
   } else {
      do.call(c, x)
   }
}

#store names of BAM fields
bam_field <- names(bam[[1]])

#go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))

#store as data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field

dim(bam_df)
#[1] 3900410      13
head(bam_df)   # 得到BAM文件


#use chr22 as an example
#how many entries on the negative strand of chr22?
table(bam_df$rname == 'chr22' & bam_df$flag == 16)
# FALSE    TRUE 
#3875997   24413

#function for checking negative strand
check_neg <- function(x){
  if (intToBits(x)[5] == 1){
    return(T)
  } else {
    return(F)
  }
}

#test neg function with subset of chr22
test <- subset(bam_df, rname == 'chr22')
dim(test)
#[1] 56426    13
table(apply(as.data.frame(test$flag), 1, check_neg))
#number same as above
#FALSE  TRUE 
#32013 24413

#function for checking positive strand
check_pos <- function(x){
  if (intToBits(x)[3] == 1){
    return(F)
  } else if (intToBits(x)[5] != 1){
    return(T)
  } else {
    return(F)
  }
}

#check pos function
table(apply(as.data.frame(test$flag), 1, check_pos))
#looks OK
#FALSE  TRUE 
#24413 32013

#store the mapped positions on the plus and minus strands
chr22_neg <- bam_df[bam_df$rname == 'chr22' &
                    apply(as.data.frame(bam_df$flag), 1, check_neg),
                    'pos'
                   ]
length(chr22_neg)
#[1] 24413
head(chr22_neg)

chr22_pos <- bam_df[bam_df$rname == 'chr22' &
                    apply(as.data.frame(bam_df$flag), 1, check_pos),
                    'pos'
                   ]
length(chr22_pos)
#[1] 32013
 
#calculate the densities
chr22_neg_density <- density(chr22_neg)
chr22_pos_density <- density(chr22_pos)
 
#display the negative strand with negative values
chr22_neg_density$y <- chr22_neg_density$y * -1
 
plot(chr22_pos_density,
     ylim = range(c(chr22_neg_density$y, chr22_pos_density$y)),
     main = "Coverage plot of mapped CAGE reads",
     xlab = "Chromosome 22",
     col = 'blue',
     lwd=2.5)
lines(chr22_neg_density, lwd=2.5, col = 'red')