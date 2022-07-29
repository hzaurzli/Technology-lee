# Goal: Explore the correlation between traits
# Input: Maize trait
# Output: Fig. 3b (v12) corrplot
# Conclusion: NA
# Author: Chia-Yi Cheng
# Last updated: 2021-01-09

# 0. Environment ---------------------------------------------------------------------------------

rm(list = setdiff(ls(), lsf.str()))
setwd("")

library(corrplot)

# 1. Load data ----------------------------------------------------------------------------------- 

raw <- read.table(file = "20142016-phenotype.txt",header=T)
head(raw)

subset=raw[,4:14] # Nutil = GrainBiomass/TotalN
head(subset)

colnames(subset)=c("Stover Biomass","Stover N","Grain Biomass","Grain N", 
                   "Total Biomass",
                   "Total N","Grain NUtE",
                   "Total Biomass Yield","Grain Yield","Stover NUtE","Total NUtE")

corr = cor(subset)

# 2. Plot --------------------------------------------------------------------------------------------

## yellow-blue
col<-colorRampPalette(c("#097054","#f7f7f7","#FFDE00"))

## Corrplot
res1 <- cor.mtest(subset, conf.level = 0.99)

corrplot(corr, method="circle",order="hclust",type = "upper",
         addrect=3,tl.cex=1.75,col=col(200),
         cl.pos="b",tl.col="grey10",tl.srt = 25, # Text label color and rotation
         addCoef.col = "grey40", # Add coefficient of correlation
         mar = c(0,0,1,0),
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