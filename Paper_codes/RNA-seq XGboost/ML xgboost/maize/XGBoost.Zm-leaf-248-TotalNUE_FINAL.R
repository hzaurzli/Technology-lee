# Goal: This script constructs XGBoost models to predict biomass (total NUE)
# Input: phenotype and gene expresion data
# Output: feature importance ranking
# Conclusion
# Author: Chia-Yi Cheng 
# Last updated: 2020-12-02

# 0. environment ----------------------------------------------------------------------------------------

rm(list = setdiff(ls(), lsf.str()))

# 1. Load data -------------------------------------------------------------------------------------

data=read.table(file="./NutriNet-Zm.TotalNUE-248.input-for-machine-learning.tsv")

## Data cleaning: Remove the genotype without replicate
data=data[grep("IRHP",rownames(data),invert=T),]

## Only keep expression as features
data=data[,-c(2:4)]

# 2. XGBoost -------------------------------------------------------------------------------------

library(xgboost)
library(data.table)
library(ggplot2)

## Setting parameters

### About the data set
n=16 # number of genotypes
c=2  # number of samples per genotype

### About the run
k=100 # number of iteration
jj=k-1

### Hyperparameters
r=40 # nrounds
colsample=0.3 #0.25
eta=0.075 #0.1 
num_parallel_tree=1
subsample=0.33 # 0.25 #1 

### About the output structure
### Need to reset the following otherwise you'll get 'out of subscript' error
y=0
obs = matrix(nrow = k*n, ncol=2)
pred = matrix(nrow = k*n, ncol=2)
rmse = vector()
train.rmse = vector()
impt.out=vector('list',k*n)

for (i in 1:n){
  for (j in 0:jj){
    
    p=c*i
    test.index <-c(p-1,p)
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
    
    #Run XGBoost
    params <- list(booster = "gbtree", 
                   objective = "reg:linear", 
                   eta=eta, 
                   gamma=150, 
                   max_depth=6, 
                   min_child_weight=1, 
                   subsample=subsample,  
                   eval_metric="rmse",
                   colsample_bytree=colsample) 
    
    set.seed(j)
    bst.val<-xgb.train( params = params, 
                        data = dtrain, 
                        nrounds = r,
                        nfold = 5, 
                        showsd = T, 
                        stratified = T, 
                        print_every_n = 10, 
                        early_stop_round = 5, 
                        watchlist = watchlist,
                        maximize = F,
                        verbose = F)
    
    y=y+1
    
    pred[y,1:c]<- predict(bst.val, dtest)
    obs[y,1:c]<-test.trait
    rmse[y]<- as.numeric(bst.val$evaluation_log[r,3])
    train.rmse[y]<-as.numeric(bst.val$evaluation_log[r,2])
    
    # extract important features
    importance_matrix <- xgb.importance(model = bst.val)
    impt.out[[y]]<-paste(importance_matrix$Feature,importance_matrix$Gain, sep = ",")
  }}

# Calculate COR for each model
pred.mat = matrix(pred, nrow=(jj+1))
obs.mat = matrix(obs, nrow=(jj+1))

COR=vector()

for (i in 1:n){
  O=c(obs.mat[,c(i,i+n)])
  P=c(pred.mat[,c(i,i+n)])
  COR[i]=cor(P,O)
}

mean(COR)   # 0.7947679 

save(data,obs,pred,rmse,train.rmse,impt.out,k,jj,n,r,c,subsample,colsample,eta,num_parallel_tree,COR,
     file="./XGBoost.Zm-leaf-248-TotalNUE-output.RData")

## Plot the output

genotype<-unlist(strsplit(rownames(data)[seq(2,by=2,32)],"_",fixed=T
))[seq(2,by=3,48)]

genotype = factor(genotype, levels =c("B73xIHP1","B73xILP1","B73xLH82","B73xMo17",
                                      "B73xMo18W","B73xOh7B","B73xPH207","B73xPHG47","B73xPHG84",
                                      "B73", "IHP1","ILP1","LH82","Mo17","PH207","PHG84"))

df <- data.frame(genotype=genotype,
                 Correlation=COR)


par(mar=c(5,6,4,2)+0.1)
par(mgp=c(5,6,4))

p<-ggplot(data=df, aes(x=genotype, y=Correlation)) +
  geom_bar(stat="identity",fill="#7fbf7b",width =0.4)+
#  labs(x="Genotype")
  theme_minimal()

p+theme(axis.text.x = element_text(angle = 90,vjust=1))

# 3. Extract importantgenes  ----------------------------------------------------------------------

library(tidyverse)
library(data.table)

load(file="XGBoost.Zm-leaf-248-TotalNUE-output.RData")

annotation<-read.csv("C:/Users/chiayi/Desktop/genome.fasta/maize/v4/B73-AGPv4.gene-symbol-description-AtHomolog.202005.tsv",
                     header = F,sep = "\t")
names(annotation)=c("Gene","Symbol","Description","Arabidopsis homolog")

s=16 # number of genotypes, n
t=100 # number of iteration, jj+1
weighted.impt=vector("list",s)

# Convert list to data frame

for (i in 1:s){
  m=(i-1)*t+1
  n=t*i
  
  d2 <- as.data.frame(do.call(rbind,flatten(impt.out[m:n]))) %>% 
    separate(.,V1, into =c("Gene","Importance"), sep=",")
  
  ## convert importance from character to numeric
  d2$Importance=as.numeric(d2$Importance)
  
  ## convert data frame to data table for easy calculation
  setDT(d2)
  d2[,sum(Importance),by=Gene]
  output=d2[,.(SUM=sum(Importance)),by=Gene][order(-SUM)]
  weighted.impt[[i]]=paste(output$Gene,output$SUM, sep = ",")
  }

# Combine the impt
d3 <- as.data.frame(do.call(rbind,flatten(weighted.impt))) %>% 
  separate(.,V1, into =c("Gene","Importance"), sep=",")

d3$Importance=as.numeric(d3$Importance)

# Calculate the sum of importance scores
setDT(d3)
d3[,sum(Importance),by=Gene]

output=d3[,.(SUM=sum(Importance)),by=Gene][order(-SUM)]
nrow(output) #248

output2=merge(output,annotation,by="Gene",all.x = T)[order(-SUM)]

# Count the frequency of each feature
output.freq=d3[,.N,by=Gene][order(-N)]

table(output.freq$N==16) #202/248
table(output.freq$N>10)  #245/248

output3=merge(output2,output.freq,by="Gene",all.x = T)[order(-SUM)]

# Thisis Supplementary Table 4 - maize
write.table(output3,row.names = F,quote = F,sep = "\t",col.names = F,
            file = "./NutriNet-Zm248logCPM-predict-TotalNUE.XGBoost-importantgene-Athomolog-frequency.tsv")

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
