# Goal: This script performs ANOVA on global phenotypes measured for all samples from 3 reps.
# Input: Arabidopsis phenotype
# Output: Fig. 2c (v12) pie chart
# Conclusion: NA
# Author: Chia-Yi Cheng
# Last updated: 2020-12-02

# 0. Environment -----------------------------------------------------------------------------------------

rm(list = setdiff(ls(), lsf.str()))
setwd("")

col=c("#66c2a4","#b2e2e2", "#238b45", "#edf8fb")

library(xtable)

# 1. Load data -------------------------------------------------------------------------------------------

phenotype<-read.table('3reps_globalphenotype_20190829.txt',header=TRUE)
colnames(phenotype)=c("Accession","N","Replica","FW (g)","DW (g)","N%","TotalN","NUE",
                      "E%","15N","NUEa","NUpE","NUpE/DW","RNAsamples","Applied N")
str(phenotype)

# 2. ANOVA ---------------------------------------------------------------------------------------------

Replica<-phenotype$Replica
Accession=phenotype$Accession
Nitrogen<-factor(phenotype$N,levels=c("L","H"))

## NUE = (DW/appliedN) 
anova_model<-aov(y~Replica+
                   Nitrogen+
                   Accession+
                   Error(Replica)+
                   Accession:Nitrogen,
                 data.frame(y=phenotype$NUEa))

summary(anova_model)
anova_model

## Generate pie chart to visualize var

var<-as.data.frame(xtable(summary(anova_model)[[2]]))
df <- data.frame(group = row.names(var),
                 value = var[,2])

df

# Use excel to plot

# 3. Session Information ----------------------------------------------------------------------------
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
#[1] xtable_1.8-4
#
#loaded via a namespace (and not attached):
#[1] compiler_4.0.3 tools_4.0.3   