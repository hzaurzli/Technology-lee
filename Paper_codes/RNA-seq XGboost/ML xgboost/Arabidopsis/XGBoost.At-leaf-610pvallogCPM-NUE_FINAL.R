# Goal: Use Arabidopsis 610 top ranked (based on the lowest p-value) N-DEGs to predict NUE
# Input: gene expression (logCPM) and trait (NUE)
# Output: model performance (correlation between actural NUE and predicted NUE) and
#         feature importance ranking
# Conclusion: the models using conserved N-DEGs outperform the ones using 610 top ranked (p-value based) N-DEGs
# Author: Chia-Yi Cheng
# Last updated: 2021-01-11

# 0. Environment ----------------------------------------------------------------------------------------

library(xgboost)
library(data.table)
library(mlr)
library(ggpubr)
library(ggplot2)

rm(list = setdiff(ls(), lsf.str()))

# 1. Load Data ------------------------------------------------------------------------------------------

data<-read.table(file = "./ML-input/NutriNet-At-veg.610pvallogCPM-NUE.input-for-machine-learning.tsv",header = T) 

## Remove unwanted features
data=data[,-c(2:4)] 

# 2. XGBoost -------------------------------------------------------------------------------------

## Setting parameters

### About the data set
n=18 # number of genotypes
c=6  # number of samples per genotype

### About the run
k=100 # number of iteration
jj=k-1

### Hyperparameters
r=50 # number of rounds
colsample=0.33
eta=0.075 #0.075 #0.1
num_parallel_tree=1
subsample=0.25

### About the output structure
## Need to reset the following otherwise you'll get 'out of subscript' error
y=0
obs = matrix(nrow = k*n, ncol=c)
pred = matrix(nrow = k*n, ncol=c)
rmse = vector()
train.rmse = vector()
impt.out=vector('list',k*n)

for (i in 1:n){
  for (j in 123:(123+jj)) {
    
    test.index <-c(i,i+18,i+36,i+54,i+72,i+90)
    testing<-data[test.index,]
    training<-data[-test.index,]
    
    #convert data frame to data table
    setDT(training) 
    setDT(testing)
    
    #using one hard encoding 
    train.trait <- training$trait 
    test.trait <- testing$trait
    new_training <- model.matrix(~.+0,data = training[,-c("trait"),with=F]) 
    new_testing <- model.matrix(~.+0,data = testing[,-c("trait"),with=F])
    
    #preparing matrix 
    dtrain <- xgb.DMatrix(data = new_training,label = train.trait) 
    dtest <- xgb.DMatrix(data = new_testing,label=test.trait)
    watchlist <- list(train=dtrain, test=dtest)
    
    #user defined evaluation metric
    cor_metric <- function(preds, dtrain) {
      labels <- getinfo(dtrain, "label")
      cor <- cor(preds,labels)
      list(metric = "cor", value = cor)}
    
    params <- list(booster = "gbtree", 
                   objective = "reg:linear", 
                   eta= eta,
                   gamma= 50,  
                   max_depth=6,
                   min_child_weight=1, 
                   eval_metric=cor_metric,
                   subsample=subsample, 
                   colsample_bytree=colsample,
                   num_parallel_tree=num_parallel_tree) 
    
    set.seed(j)
    bst.val<-xgb.train( params = params, 
                        data = dtrain, 
                        nrounds = r,
                        nfold = 5, 
                        showsd = T, 
                        stratified = T, 
                        print_every_n = 20, 
                        early_stop_round = 5, 
                        watchlist = watchlist,
                        maximize = F,
                        verbose = 0)
    
    y=y+1
    
    pred[y,1:c]<- predict(bst.val, dtest)
    obs[y,1:c]<-test.trait
    
    rmse[y]<- as.numeric(bst.val$evaluation_log[r,3])
    train.rmse[y]<-as.numeric(bst.val$evaluation_log[r,2])
    
    # extract important features
    importance_matrix <- xgb.importance(model = bst.val)
    impt.out[[y]]<-paste(importance_matrix$Feature,importance_matrix$Gain, sep = ",")
  }}

save(data,obs, pred, impt.out, rmse, train.rmse,k,n,jj,r,c,colsample,eta,num_parallel_tree,subsample,
     file="./XGBoost.At-leaf-610pvallogCPM-NUE.RData")

## Organize cor
rm(list = setdiff(ls(), lsf.str()))
setwd("C:/Users/chiayi/Desktop/NutriNet/Publication/Cheng2018/scripts/202011/Arabidopsis/")
load(file = "XGBoost.At-leaf-610pvallogCPM-NUE.RData")

pred.mat = matrix(pred, nrow=k)
obs.mat = matrix(obs, nrow=k)

# For each accession, calculate COR for each iteration 

cor.mat=matrix(nrow = jj+1,ncol = n)

for (i in 1:n){
  for (j in 1:(jj+1)){
    O=c(obs.mat[j,c(i,i+n,i+2*n,i+3*n,i+4*n,i+5*n)])
    P=c(pred.mat[j,c(i,i+n,i+2*n,i+3*n,i+4*n,i+5*n)])
    cor.mat[j,i]=cor(P,O)
  }}

cor.mat.random=cor.mat

COR=vector()
COR=colMeans(cor.mat)

mean(COR)   #0.5958295 
median(COR) #0.6972752 

# 3. Compare the results with using conserved N-DEG
load(file="./XGBoost.At-leaf-610logCPM-NUE-r50.RData")

pred.mat = matrix(pred, nrow=(jj+1))
obs.mat = matrix(obs, nrow=(jj+1))

# Calculate COR for each accession
cor.mat=matrix(nrow = jj+1,ncol = n)
for (i in 1:n){
  for (j in 1:(jj+1)){
    O=c(obs.mat[j,c(i,i+n,i+2*n,i+3*n,i+4*n,i+5*n)])
    P=c(pred.mat[j,c(i,i+n,i+2*n,i+3*n,i+4*n,i+5*n)])
    cor.mat[j,i]=cor(P,O)
  }}

colMeans(cor.mat)
apply(cor.mat,2,sd)

mean(colMeans(cor.mat)) #0.6574543
median(colMeans(cor.mat)) #0.7527112

colMeans(cor.mat)
colMeans(cor.mat.random)

wilcox.test(colMeans(cor.mat),colMeans(cor.mat.random),paired=T,alternative=c("greater"))
#V = 161, p-value = 0.000164

# 4. Session information --------------------------------------------------------------------------
sessionInfo()
#R version 4.0.3 (2020-10-10)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows >= 8 x64 (build 9200)
#
#Matrix products: default
#
#locale:
#[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
#[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
#[5] LC_TIME=English_United States.1252    
#
#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     
#
#other attached packages:
# [1] forcats_0.5.0     stringr_1.4.0     dplyr_1.0.2       purrr_0.3.4       readr_1.4.0      
# [6] tidyr_1.1.2       tibble_3.0.4      tidyverse_1.3.0   ggpubr_0.4.0      ggplot2_3.3.2    
#[11] mlr_2.18.0        ParamHelpers_1.14 data.table_1.13.2 xgboost_1.2.0.1  
#
#loaded via a namespace (and not attached):
# [1] Rcpp_1.0.5        lubridate_1.7.9   lattice_0.20-41   assertthat_0.2.1  digest_0.6.27    
# [6] R6_2.5.0          cellranger_1.1.0  backports_1.1.10  reprex_0.3.0      httr_1.4.2       
#[11] pillar_1.4.6      rlang_0.4.8       curl_4.3          readxl_1.3.1      rstudioapi_0.11  
#[16] car_3.0-10        Matrix_1.2-18     checkmate_2.0.0   labeling_0.4.2    splines_4.0.3    
#[21] foreign_0.8-80    munsell_0.5.0     broom_0.7.2       compiler_4.0.3    modelr_0.1.8     
#[26] pkgconfig_2.0.3   BBmisc_1.11       tidyselect_1.1.0  rio_0.5.16        fansi_0.4.1      
#[31] crayon_1.3.4      dbplyr_2.0.0      withr_2.3.0       grid_4.0.3        jsonlite_1.7.1   
#[36] gtable_0.3.0      lifecycle_0.2.0   DBI_1.1.0         magrittr_1.5      scales_1.1.1     
#[41] zip_2.1.1         cli_2.1.0         stringi_1.5.3     carData_3.0-4     farver_2.0.3     
#[46] ggsignif_0.6.0    fs_1.5.0          parallelMap_1.5.0 xml2_1.3.2        ellipsis_0.3.1   
#[51] generics_0.1.0    vctrs_0.3.4       openxlsx_4.2.3    fastmatch_1.1-0   tools_4.0.3      
#[56] glue_1.4.2        hms_0.5.3         abind_1.4-5       parallel_4.0.3    survival_3.2-7   
#[61] colorspace_1.4-1  rstatix_0.6.0     rvest_0.3.6       haven_2.3.1     

