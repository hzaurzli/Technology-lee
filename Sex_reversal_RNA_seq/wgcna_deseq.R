library(DESeq2)

data = read.table("/home/ug0007/rzli/fish/tmp-zebra.txt",sep = '\t',header = T,row.names = 1)
data = data[-c(34121:34125),-c(1:6,102)]
colnames(data) = c(
  paste(rep("DMSO_8_F"),1:4,sep = "_"),paste(rep("DMSO_8_M"),1:4,sep = "_"),
  paste(rep("EME2_8_F"),1:4,sep = "_"),paste(rep("EME2_8_M"),1:4,sep = "_"),
  paste(rep("EM_8_F"),1:4,sep = "_"),paste(rep("EM_8_M"),1:4,sep = "_"),
  paste(rep("E2_8_F"),1:4,sep = "_"),paste(rep("E2_8_M"),1:4,sep = "_"),
  paste(rep("DMSO_64_F"),1:4,sep = "_"),paste(rep("DMSO_64_M"),1:4,sep = "_"),
  paste(rep("EME2_64_F"),1:4,sep = "_"),paste(rep("EME2_64_M"),1:4,sep = "_"),
  paste(rep("EM_64_F"),1:4,sep = "_"),paste(rep("EM_64_M"),1:4,sep = "_"),
  paste(rep("E2_64_F"),1:4,sep = "_"),paste(rep("E2_64_M"),1:4,sep = "_"),
  paste(rep("DMSO_96_F"),1:4,sep = "_"),paste(rep("DMSO_96_M"),1:4,sep = "_"),
  paste(rep("EME2_96_F"),1:4,sep = "_"),paste(rep("EME2_96_M"),1:4,sep = "_"),
  paste(rep("EM_96_F"),1:4,sep = "_"),paste(rep("EM_96_M"),1:4,sep = "_"),
  paste(rep("E2_96_F"),1:4,sep = "_"),paste(rep("E2_96_M"),1:4,sep = "_")
)

data_64 = data[,33:64]

data_64_use = data_64[-c(32521:32525),-c(20,7,8,29,32,13,14,21,24)]
data_64_my = data_64_use[,c(1:6,13:17)]

treat_64 = factor(c(rep("DMSO",6),rep("EM",5)))
sex_64 = factor(c(rep("F",4),rep("M",2),rep("F",3),rep("M",2)))

dds_64 <- DESeqDataSetFromMatrix(data_64_my, DataFrame(treat_64,sex_64), design= ~ sex_64 + treat_64 + sex_64:treat_64)
#dds_64 <- DESeq(dds_64,parallel = T)
dds_64 <- DESeq(dds_64)
resdata_64 <- as.data.frame(counts(dds_64, normalized=TRUE))
#res_64 <- results(dds_64)
#res_64 <- results(dds_64,contrast = list(c("treat_64_EM_vs_DMSO")))
#res_64 <- results(dds_64,contrast = list(c("sex_64_M_vs_F","sex_64M.treat_64EM")))
#res_64 <- results(dds_64,contrast = list(c("sex_64_M_vs_F")))

res_64 <- results(dds_64,contrast = list(c("treat_64_EM_vs_DMSO")))

res_64 = res_64[order(res_64$padj),]
res_64 = as.data.frame(res_64)

res_64_up=res_64[res_64$padj<0.05 & res_64$log2FoldChange>1,]
res_64_down=res_64[res_64$padj<0.05 & res_64$log2FoldChange<(-1),]
res_64_up=na.omit(res_64_up)
res_64_down=na.omit(res_64_down)
res_64_up$geneID = row.names(res_64_up)
res_64_down$geneID = row.names(res_64_down)

EM_F_DMSO_F_64 = rbind(res_64_up,res_64_down)

#####################################################
treat_64 = factor(c(rep("DMSO",6),rep("EM",5)))
sex_64 = factor(c(rep("F",4),rep("M",2),rep("F",3),rep("M",2)))

sex_64 <- relevel(sex_64, "M")

dds_64 <- DESeqDataSetFromMatrix(data_64_my, DataFrame(treat_64,sex_64), design= ~ sex_64 + treat_64 + sex_64:treat_64)
dds_64 <- DESeq(dds_64)
resdata_64 <- as.data.frame(counts(dds_64, normalized=TRUE))
#res_64 <- results(dds_64)
#res_64 <- results(dds_64,contrast = list(c("treat_64_EM_vs_DMSO")))
#res_64 <- results(dds_64,contrast = list(c("sex_64_M_vs_F","sex_64M.treat_64EM")))
#res_64 <- results(dds_64,contrast = list(c("sex_64_F_vs_M")))

res_64 <- results(dds_64,contrast = list(c("treat_64_EM_vs_DMSO")))

res_64 = res_64[order(res_64$padj),]
res_64 = as.data.frame(res_64)

res_64_up=res_64[res_64$padj<0.05 & res_64$log2FoldChange>1,]
res_64_down=res_64[res_64$padj<0.05 & res_64$log2FoldChange<(-1),]
res_64_up=na.omit(res_64_up)
res_64_down=na.omit(res_64_down)
res_64_up$geneID = row.names(res_64_up)
res_64_down$geneID = row.names(res_64_down)

EM_M_DMSO_M_64 = rbind(res_64_up,res_64_down)

#####################################################################


data_32 = read.table("/home/ug0007/rzli/fish/zebrafish_32.txt",sep = '\t',header = T,row.names = 1)
data_32 = data_32[-c(32521:34125),c(1:4,9:12)]

colnames(data_32) = c(
  paste(rep("DMSO_32_F"),1:2,sep = "_"),paste(rep("DMSO_32_M"),1:2,sep = "_"),
  paste(rep("EM_32_F"),1:2,sep = "_"),paste(rep("EM_32_M"),1:2,sep = "_")
)


treat_32 = factor(c(rep("DMSO",4),rep("EM",4)))
sex_32 = factor(c(rep("F",2),rep("M",2),rep("F",2),rep("M",2)))

#sex_32 <- relevel(sex_32, "M")

dds_32 <- DESeqDataSetFromMatrix(data_32, DataFrame(treat_32,sex_32), design= ~ sex_32 + treat_32 + sex_32:treat_32)
dds_32 <- DESeq(dds_32)

#res_32 <- results(dds_32,contrast = list(c("treat_32_EM_vs_DMSO")))
#res_32 <- results(dds_32,contrast = list(c("sex_32_M_vs_F","sex_32M.treat_32EM")))
#res_32 <- results(dds_32,contrast = list(c("sex_32_M_vs_F")))

res_32 <- results(dds_32,contrast = list(c("treat_32_EM_vs_DMSO")))

res_32 = res_32[order(res_32$padj),]
res_32 = as.data.frame(res_32)

res_32_up=res_32[res_32$padj<0.05 & res_32$log2FoldChange>1,]
res_32_down=res_32[res_32$padj<0.05 & res_32$log2FoldChange<(-1),]
res_32_up=na.omit(res_32_up)
res_32_down=na.omit(res_32_down)
res_32_up$geneID = row.names(res_32_up)
res_32_down$geneID = row.names(res_32_down)

EM_F_DMSO_F_32 = rbind(res_32_up,res_32_down)

####################################################

data_32 = read.table("/home/ug0007/rzli/fish/zebrafish_32.txt",sep = '\t',header = T,row.names = 1)
data_32 = data_32[-c(32521:34125),c(1:4,9:12)]

colnames(data_32) = c(
  paste(rep("DMSO_32_F"),1:2,sep = "_"),paste(rep("DMSO_32_M"),1:2,sep = "_"),
  paste(rep("EM_32_F"),1:2,sep = "_"),paste(rep("EM_32_M"),1:2,sep = "_")
)


treat_32 = factor(c(rep("DMSO",4),rep("EM",4)))
sex_32 = factor(c(rep("F",2),rep("M",2),rep("F",2),rep("M",2)))

sex_32 <- relevel(sex_32, "M")

dds_32 <- DESeqDataSetFromMatrix(data_32, DataFrame(treat_32,sex_32), design= ~ sex_32 + treat_32 + sex_32:treat_32)
dds_32 <- DESeq(dds_32)

#res_32 <- results(dds_32,contrast = list(c("treat_32_EM_vs_DMSO")))
#res_32 <- results(dds_32,contrast = list(c("sex_32_M_vs_F","sex_32M.treat_32EM")))
#res_32 <- results(dds_32,contrast = list(c("sex_32_M_vs_F")))

res_32 <- results(dds_32,contrast = list(c("treat_32_EM_vs_DMSO")))

res_32 = res_32[order(res_32$padj),]
res_32 = as.data.frame(res_32)

res_32_up=res_32[res_32$padj<0.05 & res_32$log2FoldChange>1,]
res_32_down=res_32[res_32$padj<0.05 & res_32$log2FoldChange<(-1),]
res_32_up=na.omit(res_32_up)
res_32_down=na.omit(res_32_down)
res_32_up$geneID = row.names(res_32_up)
res_32_down$geneID = row.names(res_32_down)

EM_M_DMSO_M_32 = rbind(res_32_up,res_32_down)
####################################################
data_DMSO = read.table("/home/ug0007/rzli/fish/lav_DMSO.txt",sep = '\t',header = T,row.names = 1)
data_DMSO = data_DMSO[-c(32521:32525),-17]
data_EM = read.table("/home/ug0007/rzli/fish/lav_EM.txt",sep = '\t',header = T,row.names = 1)
data_EM = data_EM[-c(32521:32525),-18]

data_lav = cbind(data_DMSO,data_EM)
data_lav = data_lav[,c(12:16,29:33)]

day_lav = factor(c(rep("32",3),rep("64",2),rep("32",3),rep("64",2)))
treat_lav = factor(c(rep("DMSO",5),rep("EM",5)))

dds_lav <- DESeqDataSetFromMatrix(data_lav, DataFrame(treat_lav,day_lav), design= ~ day_lav + treat_lav + day_lav:treat_lav)
dds_lav <- DESeq(dds_lav)


#res_lav <- results(dds_lav)
#res_lav <- results(dds_lav,contrast = list(c("treat_lav_EM_vs_DMSO")))
#res_lav <- results(dds_lav,contrast = list(c("sex_lav_M_vs_F","sex_lavM.treat_lavEM")))
#res_lav <- results(dds_lav,contrast = list(c("sex_lav_M_vs_F")))

res_lav <- results(dds_lav,contrast = list(c("treat_lav_EM_vs_DMSO")))

res_lav = res_lav[order(res_lav$padj),]
res_lav = as.data.frame(res_lav)

res_lav_up=res_lav[res_lav$padj<0.05 & res_lav$log2FoldChange>1,]
res_lav_down=res_lav[res_lav$padj<0.05 & res_lav$log2FoldChange<(-1),]
res_lav_up=na.omit(res_lav_up)
res_lav_down=na.omit(res_lav_down)
res_lav_up$geneID = row.names(res_lav_up)
res_lav_down$geneID = row.names(res_lav_down)

lav_EM_DMSO_32 = rbind(res_lav_up,res_lav_down)

############################################
day_lav = factor(c(rep("32",3),rep("64",2),rep("32",3),rep("64",2)))
treat_lav = factor(c(rep("DMSO",5),rep("EM",5)))

day_lav <- relevel(day_lav, "64")

dds_lav <- DESeqDataSetFromMatrix(data_lav, DataFrame(treat_lav,day_lav), design= ~ day_lav + treat_lav + day_lav:treat_lav)
dds_lav <- DESeq(dds_lav)


#res_lav <- results(dds_lav)
#res_lav <- results(dds_lav,contrast = list(c("treat_lav_EM_vs_DMSO")))
#res_lav <- results(dds_lav,contrast = list(c("sex_lav_M_vs_F","sex_lavM.treat_lavEM")))
#res_lav <- results(dds_lav,contrast = list(c("sex_lav_M_vs_F")))

res_lav <- results(dds_lav,contrast = list(c("treat_lav_EM_vs_DMSO")))

res_lav = res_lav[order(res_lav$padj),]
res_lav = as.data.frame(res_lav)

res_lav_up=res_lav[res_lav$padj<0.05 & res_lav$log2FoldChange>1,]
res_lav_down=res_lav[res_lav$padj<0.05 & res_lav$log2FoldChange<(-1),]
res_lav_up=na.omit(res_lav_up)
res_lav_down=na.omit(res_lav_down)
res_lav_up$geneID = row.names(res_lav_up)
res_lav_down$geneID = row.names(res_lav_down)

lav_EM_DMSO_64 = rbind(res_lav_up,res_lav_down)


a=EM_F_DMSO_F_64[,7]
b=EM_M_DMSO_M_64[,7]
c=EM_F_DMSO_F_32[,7]
d=EM_M_DMSO_M_32[,7]
e=lav_EM_DMSO_32[,7]
f=lav_EM_DMSO_64[,7]


select_gene = union(a,b)
select_gene = union(select_gene,c)
select_gene = union(select_gene,d)
select_gene = union(select_gene,e)
select_genes = union(select_gene,f)

exp = data.frame(select_genes)

vsd_32 <- vst(dds_32, blind = FALSE)
res_32 = data.frame(assay(vsd_32))

vsd_64 <- vst(dds_64, blind = FALSE)
res_64 = data.frame(assay(vsd_64))

vsd_lav <- vst(dds_lav, blind = FALSE)
res_lav = data.frame(assay(vsd_lav))


