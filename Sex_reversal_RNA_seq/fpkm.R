aa_96 = read.table("/home/ug0007/rzli/gene_fpkm_96",header = T)
ab_96 = data.frame(aa_96$SEX0196.1.1_FRRB19H000667.1a_quant_0,
                   aa_96$SEX0196.1.2_FRRB19H000671.1a_quant_0,
                   aa_96$SEX0196.1.3_FRRB19H000991.1a_quant_0,
                   aa_96$SEX0196.1.4_FRRB19H000992.1a_quant_0,
                   aa_96$SEX0196.2.1_FRRB19H000675.1a_quant_0,
                   aa_96$SEX0196.2.2_FRRB19H000679.1a_quant_0,
                   aa_96$SEX0196.2.3_FRRB19H000999.1a_quant_0,
                   aa_96$SEX0196.2.4_FRRB19H001000.1a_quant_0,
                   aa_96$SEX0296.1.1_FRRB19H000668.1a_quant_0,
                   aa_96$SEX0296.1.2_FRRB19H000672.1a_quant_0,
                   aa_96$SEX0296.1.3_FRRB19H000993.1a_quant_0,
                   aa_96$SEX0296.1.4_FRRB19H000994.1a_quant_0,
                   aa_96$SEX0296.2.1_FRRB19H000676.1a_quant_0,
                   aa_96$SEX0296.2.2_FRRB19H000680.1a_quant_0,
                   aa_96$SEX0296.2.3_FRRB19H001001.1a_quant_0,
                   aa_96$SEX0296.2.4_FRRB19H001002.1a_quant_0,
                   aa_96$SEX0396.1.1_FRRB19H000669.1a_quant_0,
                   aa_96$SEX0396.1.2_FRRB19H000673.1a_quant_0,
                   aa_96$SEX0396.1.3_FRRB19H000995.1a_quant_0,
                   aa_96$SEX0396.1.4_FRRB19H000996.1a_quant_0,
                   aa_96$SEX0396.2.1_FRRB19H000677.1a_quant_0,
                   aa_96$SEX0396.2.2_FRRB19H000681.1a_quant_0,
                   aa_96$SEX0396.2.3_FRRB19H001003.1a_quant_0,
                   aa_96$SEX0396.2.4_FRRB19H001004.1a_quant_0,
                   aa_96$SEX0496.1.1_FRRB19H000670.1a_quant_0,
                   aa_96$SEX0496.1.2_FRRB19H000674.1a_quant_0,
                   aa_96$SEX0496.1.3_FRRB19H000997.1a_quant_0,
                   aa_96$SEX0496.1.4_FRRB19H000998.1a_quant_0,
                   aa_96$SEX0496.2.1_FRRB19H000678.1a_quant_0,
                   aa_96$SEX0496.2.2_FRRB19H000682.1a_quant_0,
                   aa_96$SEX0496.2.3_FRRB19H001005.1a_quant_0,
                   aa_96$SEX0496.2.4_FRRB19H001006.1a_quant_0)

row.names(ab_96) = aa_96$tracking_id
colnames(ab_96) = c(paste(rep("DMSO_96_F_rep"),1:4,sep = ""),paste(rep("DMSO_96_M_rep"),1:4,sep = ""),
                    paste(rep("EME2_96_F_rep"),1:4,sep = ""),paste(rep("EME2_96_M_rep"),1:4,sep = ""),
                    paste(rep("EM_96_F_rep"),1:4,sep = ""),paste(rep("EM_96_M_rep"),1:4,sep = ""),
                    paste(rep("E2_96_F_rep"),1:4,sep = ""),paste(rep("E2_96_M_rep"),1:4,sep = ""))



aa_64 = read.table("/home/ug0007/rzli/gene_fpkm_64",header = T)
ab_64 = data.frame(aa_64$SEX0164.1.1_FRRB19H000825.1a_quant_0,
                   aa_64$SEX0164.1.2_FRRB19H001045.1a_quant_0,
                   aa_64$SEX0164.1.3_FRRB19H001046.1a_quant_0,
                   aa_64$SEX0164.1.4_FRRB19H001047.1a_quant_0,
                   aa_64$SEX0164.2.1_FRRB19H000975.1a_quant_0,
                   aa_64$SEX0164.2.2_FRRB19H000976.1a_quant_0,
                   aa_64$SEX0164.2.3_FRRB19H000977.1a_quant_0,
                   aa_64$SEX0164.2.4_FRRB19H000978.1a_quant_0,
                   aa_64$SEX0264.1.1_FRRB19H000826.1a_quant_0,
                   aa_64$SEX0264.1.2_FRRB19H001048.1a_quant_0,
                   aa_64$SEX0264.1.3_FRRB19H000965.1a_quant_0,
                   aa_64$SEX0264.1.4_FRRB19H000966.1a_quant_0,
                   aa_64$SEX0264.2.1_FRRB19H000979.1a_quant_0,
                   aa_64$SEX0264.2.2_FRRB19H000980.1a_quant_0,
                   aa_64$SEX0264.2.3_FRRB19H000981.1a_quant_0,
                   aa_64$SEX0264.2.4_FRRB19H000982.1a_quant_0,
                   aa_64$SEX0364.1.1_FRRB19H000967.1a_quant_0,
                   aa_64$SEX0364.1.2_FRRB19H000968.1a_quant_0,
                   aa_64$SEX0364.1.3_FRRB19H000969.1a_quant_0,
                   aa_64$SEX0364.1.4_FRRB19H000970.1a_quant_0,
                   aa_64$SEX0364.2.1_FRRB19H000983.1a_quant_0,
                   aa_64$SEX0364.2.2_FRRB19H000984.1a_quant_0,
                   aa_64$SEX0364.2.3_FRRB19H000985.1a_quant_0,
                   aa_64$SEX0364.2.4_FRRB19H000986.1a_quant_0,
                   aa_64$SEX0464.1.1_FRRB19H000971.1a_quant_0,
                   aa_64$SEX0464.1.2_FRRB19H000972.1a_quant_0,
                   aa_64$SEX0464.1.3_FRRB19H000973.1a_quant_0,
                   aa_64$SEX0464.1.4_FRRB19H000974.1a_quant_0,
                   aa_64$SEX0464.2.1_FRRB19H000987.1a_quant_0,
                   aa_64$SEX0464.2.2_FRRB19H000988.1a_quant_0,
                   aa_64$SEX0464.2.3_FRRB19H000989.1a_quant_0,
                   aa_64$SEX0464.2.4_FRRB19H000990.1a_quant_0)

row.names(ab_64) = aa_64$tracking_id
colnames(ab_64) = c(paste(rep("DMSO_64_F_rep"),1:4,sep = ""),paste(rep("DMSO_64_M_rep"),1:4,sep = ""),
                    paste(rep("EME2_64_F_rep"),1:4,sep = ""),paste(rep("EME2_64_M_rep"),1:4,sep = ""),
                    paste(rep("EM_64_F_rep"),1:4,sep = ""),paste(rep("EM_64_M_rep"),1:4,sep = ""),
                    paste(rep("E2_64_F_rep"),1:4,sep = ""),paste(rep("E2_64_M_rep"),1:4,sep = ""))


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
                  aa_8$SEX0208.2.4_FRRB19H001036.1a_quant_0,
                  aa_8$SEX0308.1.1_FRRB19H001021.1a_quant_0,
                  aa_8$SEX0308.1.2_FRRB19H001022.1a_quant_0,
                  aa_8$SEX0308.1.3_FRRB19H001023.1a_quant_0,
                  aa_8$SEX0308.1.4_FRRB19H001024.1a_quant_0,
                  aa_8$SEX0308.2.1_FRRB19H001037.1a_quant_0,
                  aa_8$SEX0308.2.2_FRRB19H001038.1a_quant_0,
                  aa_8$SEX0308.2.3_FRRB19H001039.1a_quant_0,
                  aa_8$SEX0308.2.4_FRRB19H001040.1a_quant_0,
                  aa_8$SEX0408.1.1_FRRB19H001025.1a_quant_0,
                  aa_8$SEX0408.1.2_FRRB19H001026.1a_quant_0,
                  aa_8$SEX0408.1.3_FRRB19H001027.1a_quant_0,
                  aa_8$SEX0408.1.4_FRRB19H001028.1a_quant_0,
                  aa_8$SEX0408.2.1_FRRB19H001041.1a_quant_0,
                  aa_8$SEX0408.2.2_FRRB19H001042.1a_quant_0,
                  aa_8$SEX0408.2.3_FRRB19H001043.1a_quant_0,
                  aa_8$SEX0408.2.4_FRRB19H001044.1a_quant_0)

row.names(ab_8) = aa_8$tracking_id
colnames(ab_8) = c(paste(rep("DMSO_8_F_rep"),1:4,sep = ""),paste(rep("DMSO_8_M_rep"),1:4,sep = ""),
                   paste(rep("EME2_8_F_rep"),1:4,sep = ""),paste(rep("EME2_8_M_rep"),1:4,sep = ""),
                   paste(rep("EM_8_F_rep"),1:4,sep = ""),paste(rep("EM_8_M_rep"),1:4,sep = ""),
                   paste(rep("E2_8_F_rep"),1:4,sep = ""),paste(rep("E2_8_M_rep"),1:4,sep = ""))

###########################8
amh = ab_8['ENSDARG00000014357',]
dmrt1 = ab_8['ENSDARG00000007349',]
gsdf = ab_8['ENSDARG00000075301',]
gata4 = ab_8['ENSDARG00000098952',]
sox9a = ab_8['ENSDARG00000003293',]
sox9b = ab_8['ENSDARG00000043923',]
cyp17a1 = ab_8['ENSDARG00000033566',]
ar = ab_8['ENSDARG00000067976',]
pdgfra = ab_8['ENSDARG00000070494',]
cyp11c1 = ab_8['NSDARG00000042014',]
esr2b = ab_8['ENSDARG00000034181',]
hsd11b2 = ab_8['ENSDARG00000001975',]
nr5a1a = ab_8['ENSDARG00000103176',]
star = ab_8['ENSDARG00000006137',]
  
dm = rbind(amh,dmrt1,gsdf,gata4,sox9a,sox9b,cyp17a1,ar,pdgfra,cyp11c1,
           esr2b,hsd11b2,nr5a1a,star)
row.names(dm) = c('amh','dmrt1','gsdf','gata4','sox9a','sox9b',
                  'cyp17a1','ar','pdgfra','cyp11c1','esr2b','hsd11b2',
                  'nr5a1a','star')

#a = pheatmap(log2(dm + 1))

###########
foxL2a = ab_8['ENSDARG00000042180',]
foxL2b = ab_8['ENSDARG00000068417',]
bmp15 = ab_8['ENSDARG00000037491',]
cyp19a1a = ab_8['ENSDARG00000041348',]
gdf9 = ab_8['ENSDARG00000003229',]
vtg1 = ab_8['ENSDARG00000092233',]
vtg2 = ab_8['ENSDARG00000055809',]
vtg3 = ab_8['ENSDARG00000016448',]
vtg7 = ab_8['ENSDARG00000092419',]
figla = ab_8['ENSDARG00000087166',]
lhx8a = ab_8['ENSDARG00000002330',]
zp2.1 = ab_8['ENSDARG00000086352',]
zp3b = ab_8['ENSDARG00000039828',]
cyp11a1 = ab_8['ENSDARG00000002347',]
esr2a = ab_8['ENSDARG00000016454',]


df = rbind(foxL2a,foxL2b,bmp15,cyp19a1a,gdf9,vtg1,vtg2,vtg3,vtg7,figla,lhx8a,zp2.1,
           zp3b,cyp11a1,esr2a)
row.names(df) = c('foxL2a','foxL2b','bmp15','cyp19a1a','gdf9','vtg1',
                  'vtg2','vtg3','vtg7','figla','lhx8a','zp2.1',
                  'zp3b','cyp11a1','esr2a')
#b = pheatmap(log2(df + 1))

library(pheatmap)

dmf8 = rbind(dm,df)
#c = pheatmap(log2(dmf8 + 1))
d = pheatmap(log2(dmf8 + 1),cluster_cols = F,cluster_rows = F)

######################64
amh = ab_64['ENSDARG00000014357',]
dmrt1 = ab_64['ENSDARG00000007349',]
gsdf = ab_64['ENSDARG00000075301',]
gata4 = ab_64['ENSDARG00000098952',]
sox9a = ab_64['ENSDARG00000003293',]
sox9b = ab_64['ENSDARG00000043923',]
cyp17a1 = ab_64['ENSDARG00000033566',]
ar = ab_64['ENSDARG00000067976',]
pdgfra = ab_64['ENSDARG00000070494',]
cyp11c1 = ab_64['NSDARG00000042014',]
esr2b = ab_64['ENSDARG00000034181',]
hsd11b2 = ab_64['ENSDARG00000001975',]
nr5a1a = ab_64['ENSDARG00000103176',]
star = ab_64['ENSDARG00000006137',]

dm = rbind(amh,dmrt1,gsdf,gata4,sox9a,sox9b,cyp17a1,ar,pdgfra,cyp11c1,
           esr2b,hsd11b2,nr5a1a,star)
row.names(dm) = c('amh','dmrt1','gsdf','gata4','sox9a','sox9b',
                  'cyp17a1','ar','pdgfra','cyp11c1','esr2b','hsd11b2',
                  'nr5a1a','star')
#a = pheatmap(log2(dm + 1))

###########
foxL2a = ab_64['ENSDARG00000042180',]
foxL2b = ab_64['ENSDARG00000068417',]
bmp15 = ab_64['ENSDARG00000037491',]
cyp19a1a = ab_64['ENSDARG00000041348',]
gdf9 = ab_64['ENSDARG00000003229',]
vtg1 = ab_64['ENSDARG00000092233',]
vtg2 = ab_64['ENSDARG00000055809',]
vtg3 = ab_64['ENSDARG00000016448',]
vtg7 = ab_64['ENSDARG00000092419',]
figla = ab_64['ENSDARG00000087166',]
lhx8a = ab_64['ENSDARG00000002330',]
zp2.1 = ab_64['ENSDARG00000086352',]
zp3b = ab_64['ENSDARG00000039828',]
cyp19a1a = ab_64['ENSDARG00000041348',]
cyp11a1 = ab_64['ENSDARG00000002347',]
esr2a = ab_64['ENSDARG00000016454',]

df = rbind(foxL2a,foxL2b,bmp15,cyp19a1a,gdf9,vtg1,vtg2,vtg3,vtg7,figla,lhx8a,zp2.1,
           zp3b,cyp11a1,esr2a)
row.names(df) = c('foxL2a','foxL2b','bmp15','cyp19a1a','gdf9','vtg1',
                  'vtg2','vtg3','vtg7','figla','lhx8a','zp2.1',
                  'zp3b','cyp11a1','esr2a')
#b = pheatmap(log2(df + 1))


dmf64 = rbind(dm,df)
#c = pheatmap(log2(dmf64+ 1))
d = pheatmap(log2(dmf64 + 1),cluster_cols = F,cluster_rows = F)

#################96
amh = ab_96['ENSDARG00000014357',]
dmrt1 = ab_96['ENSDARG00000007349',]
gsdf = ab_96['ENSDARG00000075301',]
gata4 = ab_96['ENSDARG00000098952',]
sox9a = ab_96['ENSDARG00000003293',]
sox9b = ab_96['ENSDARG00000043923',]
cyp17a1 = ab_96['ENSDARG00000033566',]
ar = ab_96['ENSDARG00000067976',]
pdgfra = ab_96['ENSDARG00000070494',]
cyp11c1 = ab_96['NSDARG00000042014',]
esr2b = ab_96['ENSDARG00000034181',]
hsd11b2 = ab_96['ENSDARG00000001975',]
nr5a1a = ab_96['ENSDARG00000103176',]
star = ab_96['ENSDARG00000006137',]


dm = rbind(amh,dmrt1,gsdf,gata4,sox9a,sox9b,cyp17a1,ar,pdgfra,cyp11c1,
           esr2b,hsd11b2,nr5a1a,star)
row.names(dm) = c('amh','dmrt1','gsdf','gata4','sox9a','sox9b',
                  'cyp17a1','ar','pdgfra','cyp11c1','esr2b','hsd11b2',
                  'nr5a1a','star')

#a = pheatmap(log2(dm + 1))

###########
foxL2a = ab_96['ENSDARG00000042180',]
foxL2b = ab_96['ENSDARG00000068417',]
bmp15 = ab_96['ENSDARG00000037491',]
cyp19a1a = ab_96['ENSDARG00000041348',]
gdf9 = ab_96['ENSDARG00000003229',]
vtg1 = ab_96['ENSDARG00000092233',]
vtg2 = ab_96['ENSDARG00000055809',]
vtg3 = ab_96['ENSDARG00000016448',]
vtg7 = ab_96['ENSDARG00000092419',]
figla = ab_96['ENSDARG00000087166',]
lhx8a = ab_96['ENSDARG00000002330',]
zp2.1 = ab_96['ENSDARG00000086352',]
zp3b = ab_96['ENSDARG00000039828',]
cyp19a1a = ab_96['ENSDARG00000041348',]
cyp11a1 = ab_96['ENSDARG00000002347',]
esr2a = ab_96['ENSDARG00000016454',]

df = rbind(foxL2a,foxL2b,bmp15,cyp19a1a,gdf9,vtg1,vtg2,vtg3,vtg7,figla,lhx8a,zp2.1,
           zp3b,cyp11a1,esr2a)
row.names(df) = c('foxL2a','foxL2b','bmp15','cyp19a1a','gdf9','vtg1',
                  'vtg2','vtg3','vtg7','figla','lhx8a','zp2.1',
                  'zp3b','cyp11a1','esr2a')
#b = pheatmap(log2(df + 1))


dmf96= rbind(dm,df)
#c = pheatmap(log2(dmf96 + 1))
d = pheatmap(log2(dmf96+ 1),cluster_cols = F,cluster_rows = F)


dmf = cbind(dmf8,dmf64,dmf96)
e=pheatmap(log2(dmf+ 1),cluster_cols = F,cluster_rows = F)

###########

library(pheatmap)
library(Rtsne)
library(ggplot2)
library(ggfortify)
library(mvtnorm)
dj = dmf[,c(1:8,17:24)]
a = c(rep("DMSO_96_F",4),rep("DMSO_96_M",4),
      rep("EME2_96_F",4),rep("EME2_96_M",4),
      rep("EM_96_F",4),rep("EM_96_M",4),
      rep("E2_96_F",4),rep("E2_96_M",4))
b = c(rep("DMSO_96_F",4),rep("DMSO_96_M",4),
      rep("EM_96_F",4),rep("EM_96_M",4))

dr = rbind(as.data.frame(dj),group= b) %>%
  t()
pca = prcomp(t(dj))
autoplot(pca,data = dr,colour = 'group',size=2 )+theme_bw()


##########
rr = read.table("/home/ug0007/rzli/bio/public14.txt")

BP_up <- enrichGO(rr$V1,
                  "org.Dr.eg.db",
                  keyType="SYMBOL",
                  ont="BP"
)

barplot(BP_up, showCategory=20,title = "Biological process")
BB = BP_up@result

write.csv(BB,"/home/ug0007/rzli/bio/public.csv")


dmso_8 = ab_8[,1:8]
dmso_64 = ab_64[,1:8]
dmso_96 = ab_96[,1:8]

all = data.frame(dmso_8[,1:4],dmso_64[,1:4],dmso_96[,1:4],
                 dmso_8[,5:8],dmso_64[,5:8],dmso_96[,5:8])

pheatmap::pheatmap(cor(all,method = 'spearman'),display_numbers = T,cluster_rows = F,cluster_cols = F)


all_f = mean(cor(all[,1:12],method = 'spearman'))

all_m = mean(cor(all[,13:24],method = 'spearman'))



sex = c("female", "male","all")
data <- data.frame(mean = c(0.9135705, 0.7809227,0.7858117), sd = c(0.03025669, 0.1145665,0.1060524))
data$sex = sex


library(ggplot2)
# Fefault bar plot
p1 <- ggplot(data, aes(x = sex, y = mean)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, position = position_dodge(0.9))


all0 = all[,-c(19,21,24)]
all0_f = mean(cor(all0[,1:12],method = 'spearman'))

all0_m = mean(cor(all0[,13:21],method = 'spearman'))

sex = c("female", "male","all")
data <- data.frame(mean = c(0.9135705, 0.8556273,0.8241387), sd = c(0.03025669, 0.06763344,0.08365905))
data$sex = sex

library(ggplot2)
# Fefault bar plot
p1 <- ggplot(data, aes(x = sex, y = mean)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, position = position_dodge(0.9))

