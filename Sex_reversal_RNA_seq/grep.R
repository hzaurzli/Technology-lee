grep("ENSDARG00000041348",res_coding_up$geneID)
grep("ENSDARG00000041348",res_coding_down$geneID)

grep("ENSDARG00000092233",res_coding_up$geneID)
grep("ENSDARG00000092233",res_coding_down$geneID)

grep("ENSDARG00000055809",res_coding_up$geneID)
grep("ENSDARG00000055809",res_coding_down$geneID)

grep("ENSDARG00000016448",res_coding_up$geneID)
grep("ENSDARG00000016448",res_coding_down$geneID)

grep("ENSDARG00000078429",res_coding_up$geneID)
grep("ENSDARG00000078429",res_coding_down$geneID)

grep("ENSDARG00000092126",res_coding_up$geneID)
grep("ENSDARG00000092126",res_coding_down$geneID)

grep("ENSDARG00000016825",res_coding_up$geneID)
grep("ENSDARG00000016825",res_coding_down$geneID)

grep("ENSDARG00000092419",res_coding_up$geneID)
grep("ENSDARG00000092419",res_coding_down$geneID)

grep("ENSDARG00000014357",res_coding_up$geneID)
grep("ENSDARG00000014357",res_coding_down$geneID)






dds <- makeExampleDESeqDataSet(n=100,m=18)
dds$genotype <- factor(rep(rep(c("I","II","III"),each=3),2))
design(dds) <- ~ genotype + condition + genotype:condition
dds <- DESeq(dds)
resultsNames(dds)


gg = data.frame(results(dds, contrast= c("condition_B_vs_A","genotypeIII.conditionB") ))
nn = data.frame(results(dds, contrast=list( c("condition_B_vs_A"),c("genotypeIII.conditionB") )))
dd = data.frame(results(dds, contrast=list( c("condition_B_vs_A","genotypeIII.conditionB") )))

xx = data.frame(results(dds, name="genotypeIII.conditionB"))

data = data.frame(a = c(1,3,6,7),
                  b = c(4,1,2,3),
                  c = c(8,7,1,2),
                  d = c(5,4,3,1))
                  

p_value = apply(data, 1, function(x){
  fisher.test(matrix(c(x[1],x[2],x[3],x[4]), nrow = 2, ncol = 2))$p.value
})

data = cbind(data,p_value)


dd = data.frame(a = c(1,3,4,5),
                b = c(3,5,1,8),
                c = c(9,3,1,2))
cc = data.frame(e = c(1,5,3,2),
                f = c(7,6,3,2))
cor(dd,cc)


################################
aa_64 = read.table("/home/ug0007/rzli/gene_fpkm_64",header = T)
ab_64 = data.frame(aa_64$SEX0164.1.1_FRRB19H000825.1a_quant_0,
                   aa_64$SEX0164.1.2_FRRB19H001045.1a_quant_0,
                   aa_64$SEX0164.1.3_FRRB19H001046.1a_quant_0,
                   aa_64$SEX0164.1.4_FRRB19H001047.1a_quant_0,
                   aa_64$SEX0164.2.1_FRRB19H000975.1a_quant_0,
                   aa_64$SEX0164.2.2_FRRB19H000976.1a_quant_0,
    
                   aa_64$SEX0264.1.1_FRRB19H000826.1a_quant_0,
                   aa_64$SEX0264.1.2_FRRB19H001048.1a_quant_0,
                   aa_64$SEX0264.1.3_FRRB19H000965.1a_quant_0,
                   aa_64$SEX0264.1.4_FRRB19H000966.1a_quant_0,
       
                   aa_64$SEX0264.2.3_FRRB19H000981.1a_quant_0,
                   aa_64$SEX0264.2.4_FRRB19H000982.1a_quant_0,
                   aa_64$SEX0364.1.1_FRRB19H000967.1a_quant_0,
                   aa_64$SEX0364.1.2_FRRB19H000968.1a_quant_0,
                   aa_64$SEX0364.1.3_FRRB19H000969.1a_quant_0,
                  
                   aa_64$SEX0364.2.2_FRRB19H000984.1a_quant_0,
                   aa_64$SEX0364.2.3_FRRB19H000985.1a_quant_0,
               
                   aa_64$SEX0464.1.1_FRRB19H000971.1a_quant_0,
                   aa_64$SEX0464.1.2_FRRB19H000972.1a_quant_0,
                   aa_64$SEX0464.1.3_FRRB19H000973.1a_quant_0,
                   aa_64$SEX0464.1.4_FRRB19H000974.1a_quant_0,
      
                   aa_64$SEX0464.2.2_FRRB19H000988.1a_quant_0,
                   aa_64$SEX0464.2.3_FRRB19H000989.1a_quant_0
         )

row.names(ab_64) = aa_64$tracking_id
colnames(ab_64) = c(paste(rep("DMSO_64_F_rep"),1:4,sep = ""),paste(rep("DMSO_64_M_rep"),1:2,sep = ""),
                    paste(rep("EME2_64_F_rep"),1:4,sep = ""),paste(rep("EME2_64_M_rep"),3:4,sep = ""),
                    paste(rep("EM_64_F_rep"),1:3,sep = ""),paste(rep("EM_64_M_rep"),2:3,sep = ""),
                    paste(rep("E2_64_F_rep"),1:4,sep = ""),paste(rep("E2_64_M_rep"),2:3,sep = ""))


aa_8 = read.table("/home/ug0007/rzli/gene_fpkm_8",header = T)
ab_8 = data.frame(aa_8$SEX0108.1.1_FRRB19H001013.1a_quant_0,
                  aa_8$SEX0108.1.2_FRRB19H001014.1a_quant_0,
                  aa_8$SEX0108.1.3_FRRB19H001015.1a_quant_0,
                  aa_8$SEX0108.1.4_FRRB19H001016.1a_quant_0,
                  aa_8$SEX0108.2.1_FRRB19H001029.1a_quant_0,
                  aa_8$SEX0108.2.2_FRRB19H001030.1a_quant_0,
                  aa_8$SEX0108.2.3_FRRB19H001031.1a_quant_0,
                  aa_8$SEX0108.2.4_FRRB19H001032.1a_quant_0,
                  aa_8$SEX0208.1.1_FRRB19H001017.1a_quant_0,
                  aa_8$SEX0208.1.2_FRRB19H001018.1a_quant_0,
                  aa_8$SEX0208.1.3_FRRB19H001019.1a_quant_0,
                  aa_8$SEX0208.1.4_FRRB19H001020.1a_quant_0,
                  aa_8$SEX0208.2.1_FRRB19H001033.1a_quant_0,
                  aa_8$SEX0208.2.2_FRRB19H001034.1a_quant_0,
                  aa_8$SEX0208.2.3_FRRB19H001035.1a_quant_0,
                  
                  aa_8$SEX0308.1.3_FRRB19H001023.1a_quant_0,
                  aa_8$SEX0308.1.4_FRRB19H001024.1a_quant_0,
                  aa_8$SEX0308.2.1_FRRB19H001037.1a_quant_0,
                  aa_8$SEX0308.2.2_FRRB19H001038.1a_quant_0,
                  aa_8$SEX0308.2.3_FRRB19H001039.1a_quant_0,
                  aa_8$SEX0308.2.4_FRRB19H001040.1a_quant_0,
           
                  aa_8$SEX0408.1.2_FRRB19H001026.1a_quant_0,
                  aa_8$SEX0408.1.3_FRRB19H001027.1a_quant_0,
                  aa_8$SEX0408.1.4_FRRB19H001028.1a_quant_0,
                  
                  aa_8$SEX0408.2.2_FRRB19H001042.1a_quant_0,
                  aa_8$SEX0408.2.3_FRRB19H001043.1a_quant_0
                )

row.names(ab_8) = aa_8$tracking_id
colnames(ab_8) = c(paste(rep("DMSO_8_F_rep"),1:4,sep = ""),paste(rep("DMSO_8_M_rep"),1:4,sep = ""),
                   paste(rep("EME2_8_F_rep"),1:4,sep = ""),paste(rep("EME2_8_M_rep"),1:3,sep = ""),
                   paste(rep("EM_8_F_rep"),3:4,sep = ""),paste(rep("EM_8_M_rep"),1:4,sep = ""),
                   paste(rep("E2_8_F_rep"),2:4,sep = ""),paste(rep("E2_8_M_rep"),2:3,sep = ""))

###########################
p = pheatmap(cor(ab_8,method = "spearman"))
p = pheatmap(cor(ab_64,method = "spearman"))

###########################

diabloa = ab_64['ENSDARG00000104172',]
socs3a = ab_64['ENSDARG00000025428',]
BCL2 = ab_64['ENSDARG00000069282',]
bnip3lb = ab_64['ENSDARG00000028067',]
tnfaip8l2a = ab_64['ENSDARG00000075592',]
dkey = ab_64['ENSDARG00000104028',]
agr2 = ab_64['ENSDARG00000070480',]
casp22 = ab_64['ENSDARG00000091926',]
rassf6 = ab_64['ENSDARG00000000804',]
prnprs3 = ab_64['ENSDARG00000003705',]
socs3b = ab_64['ENSDARG00000026611',]
pycard = ab_64['ENSDARG00000040076',]
bik = ab_64['ENSDARG00000045549',]
rhoab = ab_64['ENSDARG00000094673',]
bmp8a = ab_64['ENSDARG00000035677',]
rps3 = ab_64['ENSDARG00000103007',]

diabloa = ab_8['ENSDARG00000104172',]
socs3a = ab_8['ENSDARG00000025428',]
BCL2 = ab_8['ENSDARG00000069282',]
bnip3lb = ab_8['ENSDARG00000028067',]
tnfaip8l2a = ab_8['ENSDARG00000075592',]
dkey = ab_8['ENSDARG00000104028',]
agr2 = ab_8['ENSDARG00000070480',]
casp22 = ab_8['ENSDARG00000091926',]
rassf6 = ab_8['ENSDARG00000000804',]
prnprs3 = ab_8['ENSDARG00000003705',]
socs3b = ab_8['ENSDARG00000026611',]
pycard = ab_8['ENSDARG00000040076',]
bik = ab_8['ENSDARG00000045549',]
rhoab = ab_8['ENSDARG00000094673',]
bmp8a = ab_8['ENSDARG00000035677',]
rps3 = ab_8['ENSDARG00000103007',]

jj_8 = rbind(diabloa,socs3a,BCL2,bnip3lb,tnfaip8l2a,dkey,agr2,casp22,
                  rassf6,prnprs3,socs3b,pycard,bik,rhoab,bmp8a,rps3)

row.names(jj_8) = c(diabloa,socs3a,BCL2,bnip3lb,tnfaip8l2a,dkey,agr2,casp22,
                   rassf6,prnprs3,socs3b,pycard,bik,rhoab,bmp8a,rps3)



pheatmap(log2(jj_8 + 1),cluster_rows = F,cluster_cols = F)


library(clusterProfiler)
data(geneList, package="DOSE")
de <- names(geneList)[abs(geneList) > 2]
ego <- enrichGO(de, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
library(enrichplot)
goplot(ego)

