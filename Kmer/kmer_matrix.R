#------------------------------------------------------------------------------
# Description: Kmer matrix explore
# Date: 2022.8.31
# Data: from SubPhaser test run "test_Arabidopsis.sh"
#------------------------------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = F)

library(dplyr)
path = 'tmp/kmer_count/'

# Read Data
# Method1
if(F){
case_info.vec <- c(2:13)

merge_df <- read.table(file = paste(path, as.character(1), '.fasta_15.fa', sep = ''), header = F)
names(merge_df) <- c('Kmer', 'Chr.1')
head(merge_df)


for (case_info in case_info.vec) {
    # print(case_info)
    total_df_filename <- paste(path, as.character(case_info), '.fasta_15.fa', sep = '')
    # print(total_df_filename)
    total_df <- read.table(file = total_df_filename, header = F)
    colnames(total_df) <- c("Kmer", paste("Chr.", as.character(case_info), sep = ""))
    # print(colnames(total_df))
    merge_df <- full_join(merge_df, total_df, by = "Kmer")
}
nrow(merge_df)
sum(is.na(merge_df))  # 45636399 NA

merge_df[1:5, 1:14]  # 'Kmer' is also a column
}

# Method2
dir(path)
length(dir(path))
merge_df <- data.frame(matrix(, nrow = 1, ncol = 2))
names(merge_df) <- c('Kmer', 'NA')
merge_df

for (i in dir(path)) {
    index = as.character(strsplit(i, split = "[.]")[[1]][1])
    print(index)
    total_df <- read.table(file = paste(path, i, sep = ''), header = F)
    colnames(total_df) <- c('Kmer', paste('Chr.', as.character(index), sep = ''))
    print(paste('Chr.', as.character(index)))
    merge_df <- full_join(merge_df, total_df, by = 'Kmer')
}
# sum(is.na(merge_df$`NA`)) # col `NA` exist a NA value in default, thus add a row
merge_df[1:5, 1:15]
# head(merge_df)


merge_df <- merge_df[-1, -2]
nrow(merge_df)
merge_df[1:5, 1:14]



# Filtering
# Note: This part is not that easy as I thought, 
# but it is still easy.
# Filter criteria: 
# 1) missing rate: 0.2, which means that a kmer could be saved if it was found in at least 80% of samples
# 2) min/max freq: 200 / 10000
# 3) min fold: 2
rownames(merge_df) <- merge_df[, 1]
merge_df <- merge_df[, -1]
head(merge_df)

# 1. simple
if(F){
# NA
# Read the author filtering criteria
# 1.
# head(na.omit(merge_df))
# merge_df <- na.omit(merge_df)
# dim(merge_df)
# 2.
# merge_df <- merge_df[complete.cases(merge_df), ]

# check
table(apply(merge_df, 1, sum, na.rm=TRUE) >= 200)
table(apply(merge_df, 1, sum, na.rm=TRUE) > 0)
}

# 2.complicated
missing.rate <- 0.2
nsample <- ncol(merge_df)
# new.merge_df <- data.frame(matrix(, nrow = 1, ncol = 13))
# names(new.merge_df) <- names(merge_df)
# new.merge_df


if(F){
case.vector <- rep(0, nrow(merge_df))
sum(is.na(merge_df[2,]))
# new.merge_df[nrow(new.merge_df)+1, ] <- merge_df[2,]

}
if(F){
kmermat.filter.1 <- function(mat) {
    # Count each lines missing rate and filtering
    for (i in 1:nrow(mat)) {
        line.missing.rate <- sum(is.na(mat[i, ]) / nsample)
        if (line.missing.rate > missing.rate | sum(mat[i, ], na.rm = T) < 200 | sum(mat[i, ], na.rm = T) > 10000 ) {
            case.vector[i] <- FALSE
        }
        else {
            # new.merge_df[nrow(new.merge_df)+1, ] <- merge_df[i,]
            case.vector[i] <- TRUE
        }
    }
}
# kmermat.filter.1(merge_df)
}


kmermat.filter.2 <- function(mat) {
    # Count each lines missing rate and filtering
    return(rowSums(mat, na.rm = T) < 200 | rowSums(mat, na.rm = T) > 10000)
}

# head(rowSums(merge_df[, c(1:5)], na.rm = T))
# head(rowSums(merge_df[, c(6:13)], na.rm = T))

c1 <- kmermat.filter.2(merge_df)
table(c1)

mat01 <- is.na(merge_df)
head(mat01)
vec01 <- (rowSums(mat01) / nsample)
head(vec01)
c2 <- vec01 < 0.2 

c3 <- data.frame(c1, c2)
head(c3)

c4 <- rowSums(c3) == 0
head(c4)

# Filter
new.merge_df <- merge_df[c4, ]
head(new.merge_df)
nrow(new.merge_df)


# Normalization
new.merge_df[is.na(new.merge_df)] <- 0
nrow(new.merge_df)
norml.merge_df <- scale(new.merge_df)
head(norml.merge_df)


# Calculate distance
library(tidyverse) 
library(cluster)
# install.packages("factoextra")
library(factoextra)

distance <- get_dist(norml.merge_df)
# fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

# K-means
k2 <- kmeans(norml.merge_df, centers = 2, nstart = 25)
str(k2)
# 1.
fviz_cluster(k2, data = norml.merge_df)
# 2. pairwise
norml.merge_df %>%
  as_tibble() %>%
  mutate(cluster = k2$cluster,
         state = row.names(norml.merge_df)) %>% 
  ggplot(aes(Chr.1, Chr.6, color = factor(cluster), label = state)) +
  geom_text()
table(k2$cluster)
#    1    2 
#   18 4846 

# PCA
# https://uc-r.github.io/pca