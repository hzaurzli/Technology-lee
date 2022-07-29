# Goal: This script performs ANOVA on global phenotypes measured for all samples from 3 reps.
# Input: 2014 phenotype
# Output: Fig. 3c (v12) pie chart
# Conclusion: NA
# Author: Chia-Yi Cheng
# Last updated: 2020-12-05

# 0. Environment ---------------------------------------------------------------------------------

library(xtable)

rm(list = setdiff(ls(), lsf.str()))
setwd("")

col=c("#238b45", "#66c2a4","#b2e2e2", "#edf8fb")

# 1. Load data ----------------------------------------------------------------------------------- 

raw <- read.table(file = "20142016-phenotype.txt",header=T)
raw$Year=factor(raw$Year)
pheno=raw

## Extract 2014 data #########################
pheno2014<-subset(pheno, pheno$Year=="2014")
str(pheno2014)

# 2. ANOVA --------------------------------------------------------------------------------------

pheno=pheno2014
Genotype<-pheno$Genotype
Nitrogen<-pheno$NRate
Trait<-pheno$TotalNUE

## ANOVA
anova_model<-aov(lm(y~Genotype*Nitrogen,
                    data=data.frame(y=Trait,Genotype=Genotype,Nitrogen=Nitrogen)))

summary(anova_model)
anova_model

## Generate pie chart to visualize var
var<-as.data.frame(xtable(summary(anova_model)))

df <- data.frame(
  group = row.names(var),
  value = var[,2]
)

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