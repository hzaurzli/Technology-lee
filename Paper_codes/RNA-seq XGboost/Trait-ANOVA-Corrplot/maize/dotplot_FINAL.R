# Goal: This script plots the NUE for 18 maize genotypes
# Input: NUE measurement
# Output: Figure 3A (v12) dotplot 
# Conclusion NA
# Author: Chia-Yi Cheng
# Last updated: 2021-01-09

# 0. Environment ----------------------------------------------------------------------------------------

rm(list = setdiff(ls(), lsf.str()))
setwd("")

library(ggplot2)

# 1. Load data ------------------------------------------------------------------------------------------

input<- read.table(file = "20142016-phenotype.txt",header=T)
colnames(input)[3]=c("Nitrogen")
filter<-input[grep("medN",input$Nitrogen,invert = T),] # exclude the data points with medN
phenotype=filter

## Assign factor
phenotype$Year<-as.factor(phenotype$Year)
phenotype$Nitrogen<-relevel(phenotype$Nitrogen,ref="L")
str(phenotype)

## Set up new order
new_order<-with(phenotype, reorder(Genotype , TotalNUE, mean))

# 2. Plot ------------------------------------------------------------------------------------------
## Total NUE = Total Biomass / Total N
TotalNUtE<-data.frame(phenotype[,c(1:3,14)])

ggplot(TotalNUtE, aes(x=new_order, y=TotalNUE,fill=Nitrogen)) +
  geom_jitter(
    aes(shape=Year, color = Nitrogen), 
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 1.5
  ) +
  # Either the following or geno_boxplot
  #stat_summary(
  #  aes(color = Nitrogen),
  #  fun.data=mean_se,  fun.args = list(mult=1), 
  #  geom = "pointrange",  size = 0.9,
  #  position = position_dodge(0.8)
  #)+
  # not both
  geom_boxplot(aes(color=Nitrogen), fill=NA,outlier.color = NA, 
               position = position_dodge(width=0.9))+
  scale_color_manual(values =  c("grey40","#ffc815"))+ #c("#67a9cf","#ef8a62")
  scale_shape_manual(values=c(19,10,4))+
  labs(y = "Total NUtE (Total biomass/Total N)") + # change y axis label
  labs(x="Genotype")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.35))+ # adjust the angle for tick label
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  theme(legend.position="right")+
  theme(legend.text = element_text(size=12))+
  theme(legend.title = element_text(size=12))

# 3. Session Information -------------------------------------------------------------------------------
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
#[1] ggplot2_3.3.2
#
#loaded via a namespace (and not attached):
# [1] digest_0.6.27    withr_2.3.0      dplyr_1.0.2      crayon_1.3.4     grid_4.0.3       R6_2.5.0        
# [7] lifecycle_0.2.0  gtable_0.3.0     magrittr_1.5     scales_1.1.1     pillar_1.4.6     rlang_0.4.8     
#[13] farver_2.0.3     rstudioapi_0.11  generics_0.1.0   vctrs_0.3.4      ellipsis_0.3.1   labeling_0.4.2  
#[19] tools_4.0.3      glue_1.4.2       purrr_0.3.4      munsell_0.5.0    compiler_4.0.3   pkgconfig_2.0.3 
#[25] colorspace_1.4-1 tidyselect_1.1.0 tibble_3.0.4    