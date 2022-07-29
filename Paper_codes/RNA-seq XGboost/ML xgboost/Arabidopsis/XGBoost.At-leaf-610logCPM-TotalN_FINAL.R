# Goal: Use Arabidopsis 610 N-DEGs to predict N uptake
# Input: gene expression (logCPM) and trait (N uptake)
# Output: model performance (correlation between actural N content and predicted N content) and
#         feature importance ranking
# Conclusion
# Author: Chia-Yi Cheng 
# Last updated: 2020-12-02

# 0. environment ----------------------------------------------------------------------------------------

library(xgboost)
library(data.table)
library(mlr)
library(ggpubr)
library(ggplot2)

rm(list = setdiff(ls(), lsf.str()))

# 1. Load data -----------------------------------------------------------------------------------

data<-read.table(file = "./ML-input/NutriNet-At-veg.610logCPM-TotalN.input-for-machine-learning.tsv",header = T) 

## Remove unwanted features
data=data[,-c(2:4)] 
data$trait = 10*(data$trait)

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
     file="./Arabidopsis/XGBoost.At-leaf-610logCPM-totalN-r50.RData")

## Organize cor
load(file = "XGBoost.At-leaf-610logCPM-totalN-r50.RData")

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

COR=vector()
COR=colMeans(cor.mat)

mean(COR)   # 0.6901726 

# Plot the output
apply(cor.mat,2,sd)
sd=apply(cor.mat,2,sd)
upper<-COR+sd
lower<-COR-sd

genotype<-factor(unlist(strsplit(rownames(data)[1:18],"_",fixed=T
))[seq(1,by=3,54)])

df <- data.frame(genotype=genotype,
                 Correlation=COR)

par(mar=c(5,6,4,2)+0.1)
par(mgp=c(5,6,4))

p<-ggplot(data=df, aes(x=genotype, y=Correlation)) +
  geom_bar(stat="identity",fill="#af8dc3",
           width =0.4,position=position_dodge())+
  geom_errorbar(aes(ymin=Correlation, ymax=upper), width=.2,colour="#af8dc3",
                position=position_dodge(.9)) +
  theme_minimal()

p+theme(axis.text.x = element_text(angle = 90,vjust=1))

#dev.off()

# 3. Extract importantgenes  ----------------------------------------------------------------------

library(tidyverse)
library(data.table)

rm(list = setdiff(ls(), lsf.str()))
setwd("C:/Users/chiayi/Desktop/NutriNet/Publication/Cheng2018/scripts/202011/Arabidopsis/")
load(file="XGBoost.At-leaf-610logCPM-totalN-r50.RData")

annotation<-read.csv("C:/Users/chiayi/Desktop/genome.fasta/Arabidopsis/Arabidopsis-gene.symbol-description-202005.tsv",
                     header = F,sep = "\t")
names(annotation)=c("Gene","Symbol","Description")

s=18 # number of genotypes, n
t=100 # number of iteration, jj+1
impt=vector("list",s)

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
  ## calculate impt
  impt[[i]]=paste(output$Gene,output$SUM, sep = ",")
  }

# Combine the impt
d3 <- as.data.frame(do.call(rbind,flatten(impt))) %>% 
  separate(.,V1, into =c("Gene","Importance"), sep=",")

d3$Importance=as.numeric(d3$Importance)

# Calculate composite score of importance
setDT(d3)

d3[,sum(Importance),by=Gene]
output=d3[,.(SUM=sum(Importance)),by=Gene][order(-SUM)]
output[1:10,]
nrow(output) #610

output2=merge(output,annotation,by="Gene",all.x = T)[order(-SUM)]
output2

# Count the frequency of each feature
output.freq=d3[,.N,by=Gene][order(-N)]
output.freq
table(output.freq$N==18) #607/610
table(output.freq$N>10)  #610/610

output3=merge(output2,output.freq,by="Gene",all.x = T)[order(-SUM)]
output3

write.table(output3,quote = F,sep = "\t",row.names = F, col.names = T,
            file="./NutriNet-At610logCPM-predict-TotalN.XGBoost-importantgene-frequency-r50.tsv")

hist(output2[,2])
hist(output.freq[,2])

# 4. Session information --------------------------------------------------------------------------
sessionInfo()

#R version 3.6.3 (2020-02-29)
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
# [1] forcats_0.5.0     stringr_1.4.0     purrr_0.3.4       readr_1.3.1       tidyr_1.0.2      
# [6] tibble_3.0.1      tidyverse_1.3.0   ggpubr_0.4.0      mlr_2.18.0        ParamHelpers_1.14
#[11] xgboost_1.2.0.1   xtable_1.8-4      data.table_1.12.8 gplots_3.1.0      pheatmap_1.0.12  
#[16] dplyr_1.0.2       Rmisc_1.5         plyr_1.8.6        lattice_0.20-41   ggplot2_3.3.0    
#[21] cape_2.0.2       
#
#loaded via a namespace (and not attached):
#  [1] colorspace_1.4-1        ggsignif_0.6.0          ellipsis_0.3.0          rio_0.5.16             
#  [5] evd_2.3-3               corpcor_1.6.9           fs_1.4.1                rstudioapi_0.11        
#  [9] mice_3.11.0             farver_2.0.3            lubridate_1.7.8         fansi_0.4.1            
# [13] xml2_1.3.2              codetools_0.2-16        splines_3.6.3           doParallel_1.0.16      
# [17] robustbase_0.93-6       knitr_1.30              jsonlite_1.6.1          broom_0.7.2            
# [21] dbplyr_2.0.0            shiny_1.5.0             compiler_3.6.3          httr_1.4.2             
# [25] backports_1.1.6         assertthat_0.2.1        Matrix_1.2-18           fastmap_1.0.1          
# [29] cli_2.1.0               later_1.1.0.1           htmltools_0.5.0         tools_3.6.3            
# [33] igraph_1.2.6            gtable_0.3.0            glue_1.4.0              fastmatch_1.1-0        
# [37] Rcpp_1.0.4.6            parallelMap_1.5.0       carData_3.0-4           cellranger_1.1.0       
# [41] vctrs_0.3.4             iterators_1.0.13        crosstalk_1.1.0.1       xfun_0.19              
# [45] rvest_0.3.6             openxlsx_4.2.3          mime_0.9                miniUI_0.1.1.1         
# [49] lifecycle_0.2.0         gtools_3.8.2            rstatix_0.6.0           DEoptimR_1.0-8         
# [53] MASS_7.3-51.6           scales_1.1.1            hms_0.5.3               promises_1.1.1         
# [57] parallel_3.6.3          RColorBrewer_1.1-2      BBmisc_1.11             HardyWeinberg_1.6.8    
# [61] curl_4.3                qpcR_1.4-1              regress_1.3-21          stringi_1.4.6          
# [65] foreach_1.5.1           checkmate_2.0.0         caTools_1.18.0          zip_2.1.1              
# [69] BiocParallel_1.18.1     manipulateWidget_0.10.1 truncnorm_1.0-8         shape_1.4.5            
# [73] rlang_0.4.8             pkgconfig_2.0.3         bitops_1.0-6            rgl_0.100.54           
# [77] Rsolnp_1.16             htmlwidgets_1.5.2       labeling_0.4.2          tidyselect_1.1.0       
# [81] magrittr_1.5            R6_2.5.0                generics_0.1.0          DBI_1.1.0              
# [85] pillar_1.4.6            haven_2.2.0             foreign_0.8-75          withr_2.3.0            
# [89] survival_3.1-12         abind_1.4-5             modelr_0.1.8            crayon_1.3.4           
# [93] car_3.0-10              fdrtool_1.2.15          KernSmooth_2.23-18      grid_3.6.3             
# [97] readxl_1.3.1            minpack.lm_1.2-1        reprex_0.3.0            digest_0.6.25          
#[101] webshot_0.5.2           httpuv_1.5.4            munsell_0.5.0   