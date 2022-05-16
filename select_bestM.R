#-------------------------------------------------------------------------------
# this script is used to select the best migration edge number based on the treemix results
#-------------------------------------------------------------------------------
### 安装所需要的包
# install.packages("SiZer")
# install.packages("OptM")
rm(list = ls())

### 加载库
library(OptM)
dat <- optM("4d/treemix_results_m20_bt10/")
dat <- optM("nCDS/rerun_randomseed_treemix_results_m20_bt10/")


### 判断哪一个m是最优的
plot_optM(dat, method = "Evanno", plot = TRUE, pdf = "2022-3-4-nCDS-OptM.pdf")
plot_optM(dat, method = "Evanno", plot = TRUE, pdf = "2022-3-23-nCDS-OptM.pdf")

### 绘制treemix结果图
library(RColorBrewer)
library(R.utils)
source("plotting_funcs.R")

for (i in 1:10){
  pdf(paste('4d_migration_13_bt_', i, '.pdf', sep = ""), width = 14, height = 7)
  plot_tree(paste('4d/treemix_results_m20_bt10/migration_13_bt_', i, sep=""))
  plot_resid(paste('4d/treemix_results_m20_bt10/migration_13_bt_', i, sep=""), '4d/poplist.txt')
  dev.off()
}

for (i in 1:10){
  pdf(paste('nCDS_migration_17_bt_', i, '.pdf', sep = ""), width = 14, height = 7)
  plot_tree(paste('nCDS/rerun_randomseed_treemix_results_m20_bt10/migration_17_bt_', i, sep=""))
  plot_resid(paste('nCDS/rerun_randomseed_treemix_results_m20_bt10/migration_17_bt_', i, sep=""), '4d/poplist.txt')
  dev.off()
}



### 获得所有run的方差
if(F){
  for (i in 1:30){
    get_f(paste('treemix_results_m20_bt10/migration_13_bt_', i, sep = ""))
  }
  
  for (i in 1:30){
	get_f(paste('rerun_randomseed_treemix_results_m20_bt10/migration_17_bt_', i, sep = ""))
  }
}

