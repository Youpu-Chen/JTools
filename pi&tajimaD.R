rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(rstatix)
library(ggpubr)

#-------------------------------------------------------------------------------
# load data
#-------------------------------------------------------------------------------
dat1 <- read.delim('出图 diversity and tajiamasD new/4D.NONV-SNPs.allin.CladeE_WPM.txt', sep = '\t')
dat2 <- read.delim('出图 diversity and tajiamasD new/4D.NONV-SNPs.allin.CladeW1_WPM.txt', sep = '\t')
dat4 <- read.delim('出图 diversity and tajiamasD new/4D.NONV-SNPs.allin.CladeW2_WPM.txt', sep = '\t')
fourD_sites <- rbind(dat1, dat2, dat4)


dat5 <- read.delim('出图 diversity and tajiamasD new/allSNPs.NONV-SNPs.allin.CladeE_WPM.txt', sep = '\t')
dat6 <- read.delim('出图 diversity and tajiamasD new/allSNPs.NONV-SNPs.allin.CladeW1_WPM.txt', sep = '\t')
dat8 <- read.delim('出图 diversity and tajiamasD new/allSNPs.NONV-SNPs.allin.CladeW2.WPM.txt', sep = '\t')
all_sites <- rbind(dat5, dat6, dat8)


#-------------------------------------------------------------------------------
# check
# remove the last line
#-------------------------------------------------------------------------------
unique(fourD_sites$scaff) # scaffold_{01..19}.  and there is a line called "Genome", which need to be removed
unique(fourD_sites$pop)

# create new data
fourD_sites <- fourD_sites[-which(fourD_sites$scaff == 'Genome'),] 
fourD_sites$pop <- factor(fourD_sites$pop, levels = c('CladeEallin', 'CladeW1allin', 'CladeW2'))


unique(all_sites$scaff) # scaffold_{01..19}.  and there is a line called "Genome", which need to be removed
unique(all_sites$pop)
# create new data
all_sites <- all_sites[-which(all_sites$scaff == 'Genome'),]
all_sites$pop <- factor(all_sites$pop, levels = c('CladeEall', 'CladeW1all', 'CladeW2'))


#-------------------------------------------------------------------------------
# Wilcoxon rank sum test
#-------------------------------------------------------------------------------
# ?var.test()
# ?wilcox.test()
# four-degenerate sites为dat{1,2,4}

if(F){
  # dat1 and dat2
  fourD_sites_dat1and2 <- rbind(dat1, dat2)
  fourD_sites_dat1and2$pop <- as.factor(fourD_sites_dat1and2$pop)
  fourD_sites_dat1and2 <- fourD_sites_dat1and2[which(fourD_sites_dat1and2$scaff != 'Genome'),]
  unique(fourD_sites_dat1and2$pop)
  fouD_sites_dat1and2_wilcoxtest_results <- fourD_sites_dat1and2 %>% wilcox_test(D ~ pop) %>% add_significance()
  
  # dat1 and dat4
  fourD_sites_dat1and4 <- rbind(dat1, dat4)
  fourD_sites_dat1and4$pop <- as.factor(fourD_sites_dat1and4$pop)
  fourD_sites_dat1and4 <- fourD_sites_dat1and4[which(fourD_sites_dat1and4$scaff != 'Genome'),]
  unique(fourD_sites_dat1and4$pop)
  fouD_sites_dat1and4_wilcoxtest_results <- fourD_sites_dat1and4 %>% wilcox_test(D ~ pop) %>% add_significance()
  
  # dat2 and dat4
  fourD_sites_dat2and4 <- rbind(dat2, dat4)
  fourD_sites_dat2and4$pop <- as.factor(fourD_sites_dat2and4$pop)
  fourD_sites_dat2and4 <- fourD_sites_dat2and4[which(fourD_sites_dat2and4$scaff != 'Genome'),]
  unique(fourD_sites_dat2and4$pop)
  fouD_sites_dat2and4_wilcoxtest_results <- fourD_sites_dat2and4 %>% wilcox_test(D ~ pop) %>% add_significance()
  
}
if(F){
  # var.test(dat1$D, dat2$D)
  wilcox.test(dat1$D, dat2$D, var.equal = FALSE, exact = FALSE)
  # wilcox.test(dat1$D, dat2$D, var.equal = FALSE, correct=FALSE)
  # wilcox.test(unique(dat1$D), unique(dat2$D), exact = TRUE)
  
  # var.test(dat1$D, dat4$D)
  wilcox.test(dat1$D, dat4$D, var.equal = FALSE, exact = FALSE)
  
  # var.test(dat2$D, dat4$D)
  wilcox.test(dat2$D, dat4$D, var.equal = FALSE, exact = FALSE)
}

#-------------------------------------------------------------------------------
# ggplot2
#-------------------------------------------------------------------------------
library(ggplot2)
#-----------
# 4D-Tajima
#-----------
p1 <- ggplot(data = fourD_sites, aes(x = pop, y = D)) + 
  geom_violin(lwd = 0.05) + 
  geom_boxplot(width=0.1, lwd = 0.05, outlier.size = 0.05, outlier.alpha = 0.5) + 
  theme_classic() +
  geom_signif(comparisons = list(c("CladeEallin", "CladeW1allin"),
                                 c("CladeEallin", "CladeW2"),
                                 c("CladeW1allin", "CladeW2")),
              map_signif_level=FALSE,
              textsize=1,
              test=wilcox.test,
              step_increase=0.2,
              size = 0.1) + 
  theme(axis.line = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.1),
        axis.text = element_text(size = 3),
        axis.title = element_text(size = 3))

p1
# ggsave('2022-4-27-significance-4D_TajimaD.pdf', p1, width = 8, height = 6)
# ggsave('2022-4-27-significance-4D_TajimaD.eps', p1, width = 8, height = 6)
# ggsave('2022-4-28-significance-4D_TajimaD.pdf', p1, width = 3, height = 3)
ggsave('2022-4-29-significance-4D_TajimaD.png', p1, width = 6, height = 6, dpi = 600)


#-----------
# 4D-pi
#-----------
p2 <- ggplot(data = fourD_sites, aes(x = pop, y = Diversity)) + 
  # geom_violin() + 
  geom_boxplot(width=0.5, lwd = 0.05, outlier.size = 0.05, outlier.alpha = 0.5) + 
  theme_classic() +
  geom_signif(comparisons = list(c("CladeEallin", "CladeW1allin"),
                                 c("CladeEallin", "CladeW2"),
                                 c("CladeW1allin", "CladeW2")),
              map_signif_level=FALSE,
              textsize=1,
              test=wilcox.test,
              step_increase=0.2,
              size = 0.1) + 
  theme(axis.line = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.1),
        axis.text = element_text(size = 3),
        axis.title = element_text(size = 3))


p2
# ggsave('2022-4-27-significance-4D_Diversity.pdf', p2, width = 8, height = 6)
# ggsave('2022-4-27-significance-4D_Diversity.eps', p2, width = 8, height = 6)
# ggsave('2022-4-28-significance-4D_Diversity.pdf', p2, width = 3, height = 3)
ggsave('2022-4-29-significance-4D_Diversity.png', p2, width = 6, height = 6, dpi = 600)



#-----------
# all-Tajima
#-----------
p3 <- ggplot(data = all_sites, aes(x = pop, y = D)) + 
  geom_violin(lwd = 0.05) + 
  geom_boxplot(width=0.1, lwd = 0.05, outlier.size = 0.05, outlier.alpha = 0.5) + 
  theme_classic() +
  geom_signif(comparisons = list(c("CladeEall", "CladeW1all"),
                                 c("CladeEall", "CladeW2"),
                                 c("CladeW1all", "CladeW2")),
              map_signif_level=FALSE,
              textsize=1,
              test=wilcox.test,
              step_increase=0.2,
              size = 0.1) + 
  theme(axis.line = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.1),
        axis.text = element_text(size = 3),
        axis.title = element_text(size = 3))

p3
# ggsave('2022-4-27-significance-all_TajimaD.pdf', p3, width = 8, height = 6)
# ggsave('2022-4-27-significance-all_TajimaD.eps', p3, width = 8, height = 6)
# ggsave('2022-4-28-significance-all_TajimaD.pdf', p3, width = 3, height = 3)
ggsave('2022-4-29-significance-all_TajimaD.png', p3, width = 6, height = 6, dpi = 600)


#-----------
# all-pi
#-----------
p4 <- ggplot(data = all_sites, aes(x = pop, y = Diversity)) + 
  geom_boxplot(width=0.5, lwd = 0.05, outlier.size = 0.05, outlier.alpha = 0.5) + 
  # geom_violin() + 
  theme_classic() +
  geom_signif(comparisons = list(c("CladeEall", "CladeW1all"),
                                 c("CladeEall", "CladeW2"),
                                 c("CladeW1all", "CladeW2")),
              map_signif_level=FALSE,
              textsize=1,
              test=wilcox.test,
              step_increase=0.2,
              size = 0.1) + 
  theme(axis.line = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.1),
        axis.text = element_text(size = 3),
        axis.title = element_text(size = 3))
p4
# ggsave('2022-4-27-significance-all_Diversity.pdf', p4, width = 8, height = 6)
# ggsave('2022-4-27-significance-all_Diversity.eps', p4, width = 8, height = 6)
# ggsave('2022-4-28-significance-all_Diversity.pdf', p4, width = 3, height = 3)
ggsave('2022-4-29-significance-all_Diversity.png', p4, width = 6, height = 6, dpi = 600)


#-------------------------------------------------------------------------------
# 统计对应群体的mean 和 stand deviation
#-------------------------------------------------------------------------------
fourD_summary_df_mean <- as.data.frame(fourD_sites %>% group_by(pop) %>% summarise(TajimaD_mean = mean(D),
                                                                              TajimaD_SD = sd(D),
                                                                              pi_mean = mean(Diversity),
                                                                              pi_SD = sd(Diversity)))

fourD_summary_df_median <- as.data.frame(fourD_sites %>% group_by(pop) %>% summarise(TajimaD_median = median(D),
                                                                              TajimaD_SD = sd(D),
                                                                              pi_median = median(Diversity),
                                                                              pi_SD = sd(Diversity)))

all_summary_df_mean <- as.data.frame(all_sites %>% group_by(pop) %>% summarise(TajimaD_mean = mean(D),
                                                                          TajimaD_SD = sd(D),
                                                                          pi_mean = mean(Diversity),
                                                                          pi_SD = sd(Diversity)))
all_summary_df_median <- as.data.frame(all_sites %>% group_by(pop) %>% summarise(TajimaD_median = median(D),
                                                                          TajimaD_SD = sd(D),
                                                                          pi_median = median(Diversity),
                                                                          pi_SD = sd(Diversity)))

library(openxlsx)
# sheets = list("sheet1"=fourD_summary_df_mean,
#               "sheet2"=all_summary_df_mean)
# write.xlsx(sheets, "2022-3-9-Mean&SD.xlsx")

sheets = list("sheet1"=fourD_summary_df_mean,
              "sheet2"=all_summary_df_mean)
write.xlsx(sheets, "2022-3-18-Median&SD.xlsx")
