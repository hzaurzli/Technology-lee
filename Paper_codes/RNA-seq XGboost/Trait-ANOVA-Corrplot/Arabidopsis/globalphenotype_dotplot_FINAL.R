# Goal: This script plots the NUE for 18 Arabidopsis accessions
# Input: NUE measurement
# Output: dotplot
# Conclusion NA
# Author: Chia-Yi Cheng
# Last updated: 2021-01-09

# 0. Environment ----------------------------------------------------------------------------------------

rm(list = setdiff(ls(), lsf.str()))
setwd("")

library(ggplot2)

# 1. Load data ------------------------------------------------------------------------------------------
phenotype<-read.table("3reps_globalphenotype_20190829.txt",header = T)

colnames(phenotype)=c("Accession","Nitrogen","Replica","FW","DW","N%","TotalN","NUtE","E%","15N",
                      "NUE","NUpE","NUpE/DW","RNAsamples","Applied N")
str(phenotype)

Replica=factor(phenotype$Replica,levels=c("one","two","three"))
Accession=phenotype$Accession
Nitrogen=factor(phenotype$Nitrogen,levels=c("L","H"))

## Set up new order based on the mean of NUE
new_order<-with(phenotype, reorder(Accession , NUE, mean))

# 2. Plot ---------------------------------------------------------------------------------------------
NUE<-data.frame(phenotype[,c(1:3,11)])
NUE$Nitrogen <- ordered(NUE$Nitrogen,levels=c("L","H"))
NUE$Replica <- ordered(NUE$Replica,levels=c("one","two","three"))

ggplot(NUE, aes(x=new_order, y=NUE,fill=Nitrogen)) +
  geom_jitter(
    aes(shape = Replica, color = Nitrogen), 
    position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.9),
    size = 1.2
  ) +
  # Either the following or geom_boxplot
  #  stat_summary(
  #    aes(color = Nitrogen),
  #    fun.data=mean_se,  fun.args = list(mult=1), 
  #    geom = "pointrange",  size = 0.9,
  #    position = position_dodge(0.8)
  #  )+
  # not both
  geom_boxplot(aes(color=Nitrogen), fill=NA,outlier.color = NA, 
               position = position_dodge(width=0.9))+
  scale_color_manual(values =  c("grey40","#ffc815"))+ #c("#67a9cf","#ef8a62")
  scale_shape_manual(values=c(1,19,4))+
  labs(y = "NUE (g dry weight/g applied N)") + # change y axis label
  labs(x="Accession")+
  theme_bw()+
  #  theme(axis.text.x = element_text(angle = 60))+ # adjust the angle for tick label
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
# [1] viridisLite_0.3.0 digest_0.6.27     withr_2.3.0       dplyr_1.0.2       crayon_1.3.4     
# [6] grid_4.0.3        R6_2.5.0          lifecycle_0.2.0   gtable_0.3.0      magrittr_1.5     
#[11] scales_1.1.1      pillar_1.4.6      rlang_0.4.8       farver_2.0.3      rstudioapi_0.11  
#[16] generics_0.1.0    vctrs_0.3.4       ellipsis_0.3.1    labeling_0.4.2    tools_4.0.3      
#[21] glue_1.4.2        purrr_0.3.4       munsell_0.5.0     compiler_4.0.3    pkgconfig_2.0.3  
#[26] colorspace_1.4-1  tidyselect_1.1.0  tibble_3.0.4     