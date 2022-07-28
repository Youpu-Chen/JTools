# This script is used to plot the miRNA frequency in each sample
rm(list = ls())
options(stringsAsFactors = FALSE)

# Dataclean
library(stringi)
bwt_output <- read.table('bwtoutput')
names(bwt_output) <- c('seqid', 'direction', 'microID', 'startpos', 'seq', 'quality', 'randomN')
# head(bwt_output)

bwt_output$freq <-  unlist(lapply(bwt_output$seqid, function(x){
    strsplit(x, split = "_")[[1]][2] }
))

bwt_output$microlen <-  unlist(lapply(bwt_output$seqid, function(x){
    strsplit(x, split = "_")[[1]][3] }
))

bwt_freq_df <- data.frame(bwt_output$freq, bwt_output$microlen)
names(bwt_freq_df) <- c('freq', 'microlen')

p_df <- aggregate(as.numeric(bwt_freq_df$freq), by=list(bwt_freq_df$microlen), FUN=sum)
names(p_df) <- c('len', 'freq')

# Plot
library(ggplot2)
p1 <- ggplot(data = p_df, aes(x=len, y=freq)) + geom_col()
p1
