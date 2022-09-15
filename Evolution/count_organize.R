#--------------------------------------------------------------------------------------------
# Description: Iso-seq transcript quantification based on featureCounts and co-analysis with  
#              OrthoFinder2 N0.tsv
# Date: 2022.9.13
#--------------------------------------------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)


# Load count matrix
# Note: It might take some times
raw.mat <- read.table('results/transcript/long.2and3.mat', sep = '\t', header = T)

# Load orthogroup gene list
ortho.info <- read.table('results/orthofinder2/orthologs_list.tsv', sep = '\t', header = F)
# head(ortho.info)
names(ortho.info) <- c('orthogroup', 'Geneid')


# Append corresponding orthogroup info to count matrix
# Note: some gene ids from GTF might not be assigned to orthogroup
merge_df <- merge(raw.mat, ortho.info, by = 'Geneid', all = T)
head(merge_df)

# Extract info
extract_df <- merge_df[, c(1, 7, 8)]
head(extract_df)
colnames(extract_df) <- c('Geneid', 'Count', 'Orthogroup')
nrow(extract_df)


# Group df, and summarize the freq
# 1. Count total read counts of each Orthogroup
if(F){
totalcount.df <- extract_df %>%
  group_by(Orthogroup) %>%
  summarize(count = sum(Count, na.rm = TRUE))
totalcount.df <- totalcount.df[which(totalcount.df$Orthogroup != "NA"), ]   # remove NA
head(totalcount.df)
}

# 2. Sort df using Orthogroup ID which'll group the same Orthogroup together
sorted.df <- extract_df[order(extract_df$Orthogroup), ]
sorted.df <- sorted.df[which(sorted.df$Orthogroup != "NA"), ]  # filter out gene which is not assigned to any Orthogroup
nrow(sorted.df)


# main function of this code
group.count <- function(line, count.vector){

    if (substr(line[1], 1, 1) == "A") {
        if (is.na(count.vector[1])) {
            count.vector[1] <- 0 + line[2]
        } else {
            count.vector[1] <- count.vector[1] + line[2]
        }
    } else if (substr(line[1], 1, 1) == "B") {
        if (is.na(count.vector[2])) {
            count.vector[2] <- 0 + line[2] 
        } else {
            count.vector[2] <- count.vector[2] + line[2]
        }
    } else if (substr(line[1], 1, 1) == "C") {
        if (is.na(count.vector[3])) {
            count.vector[3] <- 0 + line[2]
        } else {
            count.vector[3] <- count.vector[3] + line[2]
        }
    }

    return(count.vector)
}

ratio.count <- function(df){
    ortho <-  ''
    count.vector <- c(NA, NA, NA)
    row.df <- data.frame(matrix(ncol = 4, nrow = 0))
    names(row.df) <- c('Orthogroup', 'A', 'B', 'C')
    
    for (i in 1:nrow(df)) {
        # when i = 1, prog should only do the give thing
        if (i == 1) {
            print(i)
            ortho = df[i, 3]
            count.vector <- group.count(df[i, ], count.vector)

            } else if (ortho != df[i, 3]) {
                print(i)
                row.df <- rbind(row.df, c(ortho, count.vector))
                ortho = df[i, 3]
                count.vector <- c(NA, NA, NA)
                count.vector <- group.count(df[i, ], count.vector)
            
        } else {
            count.vector <- group.count(df[i, ], count.vector)
        }
    } 
    return(row.df)
}

# Run group summary
row.df <- ratio.count(sorted.df)
head(row.df)
names(row.df) <- c("Orthogroup", "A", "B", "C")
nrow(row.df)
write.csv(row.df, 'results/transcript/Hierarchical_Orthogroups.Count.csv', quote = F, row.names = F)


# Count the freq of each Ortho-alignment type
# Using python