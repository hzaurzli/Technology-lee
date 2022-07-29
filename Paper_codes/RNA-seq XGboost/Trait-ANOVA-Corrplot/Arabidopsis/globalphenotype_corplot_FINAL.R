# Goal: Explore the correlation between traits
# Input: Arabidopsis trait
# Output: Fig. 2b (v12) corrplot
# Conclusion: NA
# Author: Chia-Yi Cheng
# Last updated: 2021-01-09

# 0. Environment ---------------------------------------------------------------------------------

rm(list = setdiff(ls(), lsf.str()))
setwd("")

library(corrplot)

# 1. Load data ----------------------------------------------------------------------------------- 

phenotype<-read.table('3reps_globalphenotype_20190829.txt',header=TRUE)
colnames(phenotype)=c("Accession","N","Replica","FW (g)","DW (g)","N%","TotalN","NUE",
                      "E%","15N","NUEa","NUpE","NUpE/DW","RNAsamples","Applied N")
str(phenotype)

subset1=phenotype[,c(6,5,7,11,9,12,8)]
colnames(subset1)=c("N%","Dry weight","N uptake","NUE","E%","NUpE","NUtE")

corr = cor(subset1)
res1 <- cor.mtest(subset1, conf.level = 0.99, method="pearson")

# 2. Plot -------------------------------------------------------------------------------------------

col<-colorRampPalette(c("#097054","#f7f7f7","#FFDE00"))

corrplot(corr, method="circle",order="hclust",type = "upper",
         addrect=3,tl.cex=1.5,col=col(200),
         cl.pos="b",tl.col="grey10",tl.srt = 25, # Text label color and rotation
         addCoef.col = "grey40", # Add coefficient of correlation
         mar = c(0,0,3,0),
         insig = "blank", pch.col = "grey50",p.mat = res1[[1]],
         # hide (F) or not hide (T) the correlation coefficient on the principal diagonal
         diag = F)

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
#[1] corrplot_0.84
#
#loaded via a namespace (and not attached):
#[1] compiler_4.0.3 tools_4.0.3   