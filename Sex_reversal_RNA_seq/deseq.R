###############################
gene_lnc_id = read.table("rzli/fina.lnc.geneid",sep = ' ',header = F)
gene_lnc_id$V1 = as.character(gene_lnc_id$V1)
gene_lnc_id$V2 = as.character(gene_lnc_id$V2)
gene_lnc_id$V2 = substr(gene_lnc_id$V2,1,18)

gene_lnc_id_unique = data.frame(unique(gene_lnc_id))
gene_lnc_id_unique = data.frame(gene_lnc_id_unique[,-1])
colnames(gene_lnc_id_unique) = c("geneID")
#gene_lnc_id_unique$V2 = gsub(";","",gene_lnc_id_unique$V2)
#gene_lnc_id_unique_p = gene_lnc_id_unique[which(substr(gene_lnc_id_unique[,2],1,1) == "P"),]
#gene_lnc_id_unique_p = data.frame(gene_lnc_id_unique_p[,-1])
#colnames(gene_lnc_id_unique_p) = c("geneID")
#gene_lnc_id_unique_p$geneID = as.character(gene_lnc_id_unique_p$geneID)
#gene_lnc_id_unique_p = merge(gene_lnc_id_unique_p,id_change,by = "geneID")
#gene_lnc_id_unique = data.frame(gene_lnc_id_unique[,-1])
#colnames(gene_lnc_id_unique) = c("geneID")


gene_coding_id = read.table("rzli/fina.coding.geneid",sep = ' ',header = F)
gene_coding_id$V1 = as.character(gene_coding_id$V1)
gene_coding_id$V2 = as.character(gene_coding_id$V2)
gene_coding_id$V2 = substr(gene_coding_id$V2,1,18)

gene_coding_id_unique = data.frame(unique(gene_coding_id))
gene_coding_id_unique = data.frame(gene_coding_id_unique[,-1])
colnames(gene_coding_id_unique) = c("geneID")
#gene_coding_id_unique$V2 = gsub(";","",gene_coding_id_unique$V2)
#gene_coding_id_unique_p = gene_coding_id_unique[which(substr(gene_coding_id_unique[,2],1,1) == "P"),]
#gene_coding_id_unique_p = data.frame(gene_coding_id_unique_p[,-1])
#colnames(gene_coding_id_unique_p) = c("geneID")

#gene_coding_id_unique_e = gene_coding_id_unique[which(substr(gene_coding_id_unique[,2],1,1) == "E"),]
#gene_coding_id_unique_e = data.frame(gene_coding_id_unique_e[,-1])
#colnames(gene_coding_id_unique_e) = c("geneID")

#gene_coding_id_unique_p$geneID = as.character(gene_coding_id_unique_p$geneID)
#gene_coding_id_unique_p = merge(gene_coding_id_unique_p,id_change,by = "geneID")
#gene_coding_id_unique_p = data.frame(gene_coding_id_unique_p[,-1])
#colnames(gene_coding_id_unique_p) = c("geneID")

#gene_coding_id_unique = rbind(gene_coding_id_unique_e,gene_coding_id_unique_p)
#gene_coding_id_unique = unique(gene_coding_id_unique)

############
id_change = read.table("rzli/ensembl_pb_gene_id",sep = '\t')
id_change$V1 = as.character(id_change$V1)
id_change$V2 = as.character(id_change$V2)

colnames(id_change) = c('en_id','geneID')

####################################################
library(DESeq2)

data = read.table("rzli/tmp.txt",sep = '\t',header = T,row.names = 1)
data = data[-c(34120:34124),-53]
colnames(data) = c(paste(rep("Null_0_F"),1:2,sep = "_"),paste(rep("Null_0_M"),1:2,sep = "_"),
                      paste(rep("DMSO_8_F"),1:2,sep = "_"),paste(rep("DMSO_8_M"),1:2,sep = "_"),
                      paste(rep("EME2_8_F"),1:2,sep = "_"),paste(rep("EME2_8_M"),1:2,sep = "_"),
                      paste(rep("EM_8_F"),1:2,sep = "_"),paste(rep("EM_8_M"),1:2,sep = "_"),
                      paste(rep("E2_8_F"),1:2,sep = "_"),paste(rep("E2_8_M"),1:2,sep = "_"),
                      paste(rep("DMSO_64_F"),1:2,sep = "_"),paste(rep("DMSO_64_M"),1:2,sep = "_"),
                      paste(rep("EME2_64_F"),1:2,sep = "_"),paste(rep("EME2_64_M"),1:2,sep = "_"),
                      paste(rep("EM_64_F"),1:2,sep = "_"),paste(rep("EM_64_M"),1:2,sep = "_"),
                      paste(rep("E2_64_F"),1:2,sep = "_"),paste(rep("E2_64_M"),1:2,sep = "_"),
                      paste(rep("DMSO_96_F"),1:2,sep = "_"),paste(rep("DMSO_96_M"),1:2,sep = "_"),
                      paste(rep("EME2_96_F"),1:2,sep = "_"),paste(rep("EME2_96_M"),1:2,sep = "_"),
                      paste(rep("EM_96_F"),1:2,sep = "_"),paste(rep("EM_96_M"),1:2,sep = "_"),
                      paste(rep("E2_96_F"),1:2,sep = "_"),paste(rep("E2_96_M"),1:2,sep = "_")
                    )

nol = data.frame(rownames(data_08))
nn = cbind(data_08,nol)
colnames(nn) = c(paste(rep("DMSO_F"),1:4,sep = "_"),paste(rep("DMSO_M"),1:4,sep = "_"),
                      paste(rep("EME2_F"),1:4,sep = "_"),paste(rep("EME2_M"),1:4,sep = "_"),
                      paste(rep("EM_F"),1:4,sep = "_"),paste(rep("EM_M"),1:4,sep = "_"),
                      paste(rep("E2_F"),1:4,sep = "_"),paste(rep("E2_M"),1:4,sep = "_"),"geneID")



nol_1 = data.frame(nn[14671:34119,])
#colnames(nol_1) = c("geneID")
nol_2 = data.frame(nn[1:14670,])
#colnames(nol_2) = c("geneID")
aa = merge(nol_1,id_change,by = 'geneID')
aa = aa[,-1]
colnames(aa) = c(paste(rep("DMSO_F"),1:4,sep = "_"),paste(rep("DMSO_M"),1:4,sep = "_"),
                 paste(rep("EME2_F"),1:4,sep = "_"),paste(rep("EME2_M"),1:4,sep = "_"),
                 paste(rep("EM_F"),1:4,sep = "_"),paste(rep("EM_M"),1:4,sep = "_"),
                 paste(rep("E2_F"),1:4,sep = "_"),paste(rep("E2_M"),1:4,sep = "_"),"geneID")

id = rbind(nol_2,aa)
id = id[!duplicated(id$geneID), ]
rownames(id) = NULL
rownames(id) = id[,33]
data_08 = id[,-33]

data_08_EMF_DMSOF = data_08[,c(1:4,17:20)]

treat = factor(c("DMSO","DMSO","DMSO","DMSO","EM","EM","EM","EM"))
dds <- DESeqDataSetFromMatrix(data_08_EMF_DMSOF, DataFrame(treat), design= ~ treat)
dds <- DESeq(dds)

resdata <- as.data.frame(counts(dds, normalized=TRUE))
res <- results(dds)
res = res[order(res$padj),]
res = as.data.frame(res)
res$geneID = rownames(res)
res$geneID = as.character(res$geneID)

res_pb = res[which(substr(res[,7],1,1) == "P"),]
res_en = res[which(substr(res[,7],1,1) == "E"),]
res_pb_e = merge(res_pb,id_change,by='geneID')
res_pb_e = res_pb_e[,-1]
colnames(res_pb_e)= c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","geneID")
res_new = rbind(res_en,res_pb_e)
res_new = na.omit(res_new)
############
res_new_up=res_new[res_new$padj<0.05 & res_new$log2FoldChange>1,]
res_new_down=res_new[res_new$padj<0.05 & res_new$log2FoldChange<(-1),]
res_new_up=na.omit(res_new_up)
res_new_down=na.omit(res_new_down)

res_lnc_up = merge(gene_lnc_id_unique,res_new_up,by="geneID")
res_coding_up = merge(gene_coding_id_unique,res_new_up,by="geneID")

res_lnc_down = merge(gene_lnc_id_unique,res_new_down,by="geneID")
res_coding_down = merge(gene_coding_id_unique,res_new_down,by="geneID")


##################################3

grep("ENSDARG00000117778",res_em_dmso_m_f_lnc_up$geneID)
grep("ENSDARG00000097580",res_em_dmso_m_f_lnc_up$geneID)
grep("ENSDARG00000117226",res_em_dmso_m_f_lnc_up$geneID)
grep("ENSDARG00000117446",res_em_dmso_m_f_lnc_up$geneID)
grep("ENSDARG00000117298",res_em_dmso_m_f_lnc_up$geneID)


grep("ENSDARG00000117778",res_em_dmso_m_f_lnc_down$geneID)
grep("ENSDARG00000097580",res_em_dmso_m_f_lnc_down$geneID)
grep("ENSDARG00000117226",res_em_dmso_m_f_lnc_down$geneID)
grep("ENSDARG00000117446",res_em_dmso_m_f_lnc_down$geneID)
grep("ENSDARG00000117298",res_em_dmso_m_f_lnc_down$geneID)

grep("ENSDARG00000117778",res_coding_up$geneID)
grep("ENSDARG00000097580",res_coding_up$geneID)
grep("ENSDARG00000117226",res_coding_up$geneID)
grep("ENSDARG00000117446",res_coding_up$geneID)
grep("ENSDARG00000117298",res_coding_up$geneID)


grep("ENSDARG00000117778",res_coding_down$geneID)
grep("ENSDARG00000097580",res_coding_down$geneID)
grep("ENSDARG00000117226",res_coding_down$geneID)
grep("ENSDARG00000117446",res_coding_down$geneID)
grep("ENSDARG00000117298",res_coding_down$geneID)

#####################
#####################
#####################
library(DESeq2)

data_08 = read.table("rzli/tmp.txt",sep = '\t',header = T,row.names = 1)
data_08 = data_08[-c(34120:34124),-c(9,18,35)]
colnames(data_08) = c(paste(rep("DMSO_F"),1:4,sep = "_"),paste(rep("DMSO_M"),1:4,sep = "_"),
                      paste(rep("EME2_F"),1:4,sep = "_"),paste(rep("EME2_M"),1:4,sep = "_"),
                      paste(rep("EM_F"),1:4,sep = "_"),paste(rep("EM_M"),1:4,sep = "_"),
                      paste(rep("E2_F"),1:4,sep = "_"),paste(rep("E2_M"),1:4,sep = "_"))

nol = data.frame(rownames(data_08))
nn = cbind(data_08,nol)
colnames(nn) = c(paste(rep("DMSO_F"),1:4,sep = "_"),paste(rep("DMSO_M"),1:4,sep = "_"),
                 paste(rep("EME2_F"),1:4,sep = "_"),paste(rep("EME2_M"),1:4,sep = "_"),
                 paste(rep("EM_F"),1:4,sep = "_"),paste(rep("EM_M"),1:4,sep = "_"),
                 paste(rep("E2_F"),1:4,sep = "_"),paste(rep("E2_M"),1:4,sep = "_"),"geneID")



nol_1 = data.frame(nn[14671:34119,])
#colnames(nol_1) = c("geneID")
nol_2 = data.frame(nn[1:14670,])
#colnames(nol_2) = c("geneID")
aa = merge(nol_1,id_change,by = 'geneID')
aa = aa[,-1]
colnames(aa) = c(paste(rep("DMSO_F"),1:4,sep = "_"),paste(rep("DMSO_M"),1:4,sep = "_"),
                 paste(rep("EME2_F"),1:4,sep = "_"),paste(rep("EME2_M"),1:4,sep = "_"),
                 paste(rep("EM_F"),1:4,sep = "_"),paste(rep("EM_M"),1:4,sep = "_"),
                 paste(rep("E2_F"),1:4,sep = "_"),paste(rep("E2_M"),1:4,sep = "_"),"geneID")

id = rbind(nol_2,aa)
id = id[!duplicated(id$geneID), ]
rownames(id) = NULL
rownames(id) = id[,33]
data_08 = id[,-33]
da = data.frame(row.names(data_08))

treat = factor(c(rep("DMSO",8),rep("EME2",8),rep("EM",8),rep("E2",8)))
sex = factor(c(rep("F",4),rep("M",4),rep("F",4),rep("M",4),rep("F",4),rep("M",4),rep("F",4),rep("M",4)))
dds <- DESeqDataSetFromMatrix(data_08, DataFrame(treat,sex), design= ~ sex + treat + sex:treat)
dds <- DESeq(dds)
resdata <- as.data.frame(counts(dds, normalized=TRUE))
res <- results(dds)
res = res[order(res$padj),]
res = as.data.frame(res)

#res_em_dmso_m_f = results(dds,contrast = list(c("treat_EM_vs_DMSO","sexM.treatEM")))

#res_em_dmso_m_f = results(dds,contrast = c("treat","EM","DMSO"))
res_em_dmso_m_f = results(dds,name = c("sex_M_vs_F"))

res_em_dmso_m_f = as.data.frame(res_em_dmso_m_f)
res_em_dmso_m_f$geneID = rownames(res_em_dmso_m_f)
res_em_dmso_m_f$geneID = as.character(res_em_dmso_m_f$geneID)

#res_em_dmso_m_f_pb = res_em_dmso_m_f[which(substr(res_em_dmso_m_f[,7],1,1) == "P"),]
#res_em_dmso_m_f_en = res_em_dmso_m_f[which(substr(res_em_dmso_m_f[,7],1,1) == "E"),]
#res_em_dmso_m_f_pb_e = merge(res_em_dmso_m_f_pb,id_change,by='geneID')
#res_em_dmso_m_f_pb_e = res_em_dmso_m_f_pb_e[,-1]
#colnames(res_em_dmso_m_f_pb_e)= c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","geneID")
#res_em_dmso_m_f_new = rbind(res_em_dmso_m_f_en,res_em_dmso_m_f_pb_e)
#res_em_dmso_m_f_new = na.omit(res_em_dmso_m_f_new)

res_em_dmso_m_f_up=res_em_dmso_m_f[res_em_dmso_m_f$padj<0.05 & res_em_dmso_m_f$log2FoldChange>1,]
res_em_dmso_m_f_down=res_em_dmso_m_f[res_em_dmso_m_f$padj<0.05 & res_em_dmso_m_f$log2FoldChange<(-1),]
res_em_dmso_m_f_up=na.omit(res_em_dmso_m_f_up)
res_em_dmso_m_f_down=na.omit(res_em_dmso_m_f_down)

res_em_dmso_m_f_lnc_up = merge(gene_lnc_id_unique,res_em_dmso_m_f_up,by="geneID")
res_em_dmso_m_f_coding_up = merge(gene_coding_id_unique,res_em_dmso_m_f_up,by="geneID")

res_em_dmso_m_f_lnc_down = merge(gene_lnc_id_unique,res_em_dmso_m_f_down,by="geneID")
res_em_dmso_m_f_coding_down = merge(gene_coding_id_unique,res_em_dmso_m_f_down,by="geneID")

t_gid = res_em_dmso_m_f_up$geneID
c_gid = as.character(res_em_dmso_m_f_coding_up$geneID)
l_gid = as.character(res_em_dmso_m_f_lnc_up$geneID)
c_l_gid = c(c_gid,l_gid)
what = t_gid[which(t_gid %in% c_l_gid == FALSE)]

a1 = data.frame(what %in% gene_lnc_id_unique)
a2 = data.frame(what %in% gene_coding_id_unique)
a3 = data.frame(what %in% id_change$en_id)

library(DOSE)
library(org.Dr.eg.db)
#library(topGO)
library(clusterProfiler)
#library(pathview)
library(AnnotationHub)

BP_up <- enrichGO(res_em_dmso_m_f_new_up$geneID,
                  "org.Dr.eg.db",
                  keyType="ENSEMBL",
                  ont="BP"
)

barplot(BP_up, showCategory=20,title = "Biological process")

BP_down <- enrichGO(res_em_dmso_m_f_new_down$geneID,
                    "org.Dr.eg.db",
                    keyType="ENSEMBL",
                    ont="BP"
)

barplot(BP_down, showCategory=20,title = "Biological process")



a = pheatmap(cor(resdata,method = "spearman"),display_numbers = TRUE,number_color = "black")


##############
library(DOSE)
library(org.Dr.eg.db)
#library(topGO)
library(clusterProfiler)
#library(pathview)
library(AnnotationHub)

BP_up <- enrichGO(res_new_up$geneID,
                          "org.Dr.eg.db",
                          keyType="ENSEMBL",
                          ont="BP"
)

barplot(BP_up, showCategory=20,title = "Biological process")

BP_down <- enrichGO(res_new_down$geneID,
                  "org.Dr.eg.db",
                  keyType="ENSEMBL",
                  ont="BP"
)

barplot(BP_down, showCategory=20,title = "Biological process")

grep("ENSDARG00000041348",res_em_dmso_m_f_coding_up$geneID)
grep("ENSDARG00000041348",res_em_dmso_m_f_coding_down$geneID)

grep("ENSDARG00000092233",res_em_dmso_m_f_coding_up$geneID)
grep("ENSDARG00000092233",res_em_dmso_m_f_coding_down$geneID)

grep("ENSDARG00000055809",res_em_dmso_m_f_coding_up$geneID)
grep("ENSDARG00000055809",res_em_dmso_m_f_coding_down$geneID)

grep("ENSDARG00000016448",res_em_dmso_m_f_coding_up$geneID)
grep("ENSDARG00000016448",res_em_dmso_m_f_coding_down$geneID)

grep("ENSDARG00000078429",res_em_dmso_m_f_coding_up$geneID)
grep("ENSDARG00000078429",res_em_dmso_m_f_coding_down$geneID)

grep("ENSDARG00000092126",res_em_dmso_m_f_coding_up$geneID)
grep("ENSDARG00000092126",res_em_dmso_m_f_coding_down$geneID)

grep("ENSDARG00000016825",res_em_dmso_m_f_coding_up$geneID)
grep("ENSDARG00000016825",res_em_dmso_m_f_coding_down$geneID)

grep("ENSDARG00000092419",res_em_dmso_m_f_coding_up$geneID)
grep("ENSDARG00000092419",res_em_dmso_m_f_coding_down$geneID)

grep("ENSDARG00000014357",res_em_dmso_m_f_coding_up$geneID)
grep("ENSDARG00000014357",res_em_dmso_m_f_coding_down$geneID)


library(DESeq2)

data_08 = read.table("rzli/tmp.txt",sep = '\t',header = T,row.names = 1)
data_08 = data_08[,-c(9,18,35)]
colnames(data_08) = c(paste(rep("DMSO_F"),1:4,sep = "_"),paste(rep("DMSO_M"),1:4,sep = "_"),
                      paste(rep("EME2_F"),1:4,sep = "_"),paste(rep("EME2_M"),1:4,sep = "_"),
                      paste(rep("EM_F"),1:4,sep = "_"),paste(rep("EM_M"),1:4,sep = "_"),
                      paste(rep("E2_F"),1:4,sep = "_"),paste(rep("E2_M"),1:4,sep = "_"))



data_08_EMF_DMSOF = data_08[,c(21:24,17:20)]

treat = factor(c("DMSO","DMSO","DMSO","DMSO","EM","EM","EM","EM"))
dds <- DESeqDataSetFromMatrix(data_08_EMF_DMSOF, DataFrame(treat), design= ~ treat)
dds <- DESeq(dds)

resdata <- as.data.frame(counts(dds, normalized=TRUE))
res <- results(dds)
res = res[order(res$padj),]
res = as.data.frame(res)
res$geneID = rownames(res)
res$geneID = as.character(res$geneID)

res_pb = res[which(substr(res[,7],1,1) == "P"),]
res_en = res[which(substr(res[,7],1,1) == "E"),]
res_pb_e = merge(res_pb,id_change,by='geneID')
res_pb_e = res_pb_e[,-1]
colnames(res_pb_e)= c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","geneID")
res_new = rbind(res_en,res_pb_e)
############
res_new_up=res_new[res_new$padj<0.05 & res_new$log2FoldChange>1,]
res_new_down=res_new[res_new$padj<0.05 & res_new$log2FoldChange<(-1),]
res_new_up=na.omit(res_new_up)
res_new_down=na.omit(res_new_down)

res_lnc_up = merge(gene_lnc_id_unique,res_new_up,by="geneID")
res_coding_up = merge(gene_coding_id_unique,res_new_up,by="geneID")

res_lnc_down = merge(gene_lnc_id_unique,res_new_down,by="geneID")
res_coding_down = merge(gene_coding_id_unique,res_new_down,by="geneID")


##############
library(DOSE)
library(org.Dr.eg.db)
#library(topGO)
library(clusterProfiler)
#library(pathview)
library(AnnotationHub)

BP_up <- enrichGO(res_new_up$geneID,
                  "org.Dr.eg.db",
                  keyType="ENSEMBL",
                  ont="BP"
)

barplot(BP_up, showCategory=20,title = "Biological process")

BP_down <- enrichGO(res_new_down$geneID,
                    "org.Dr.eg.db",
                    keyType="ENSEMBL",
                    ont="BP"
)

barplot(BP_down, showCategory=20,title = "Biological process")


####################
####################
gene_lnc_id = read.table("rzli/fina.lnc.geneid",sep = ' ',header = F)
gene_lnc_id$V1 = as.character(gene_lnc_id$V1)
gene_lnc_id$V2 = as.character(gene_lnc_id$V2)
gene_lnc_id$V2 = substr(gene_lnc_id$V2,1,18)

gene_lnc_id$V2 = gsub(";","",gene_lnc_id$V2)
gene_lnc_id_unique = data.frame(unique(gene_lnc_id))
gene_lnc_id_unique = data.frame(gene_lnc_id_unique[,-1])
colnames(gene_lnc_id_unique) = c("geneID")

gene_coding_id = read.table("rzli/fina.coding.geneid",sep = ' ',header = F)
gene_coding_id$V1 = as.character(gene_coding_id$V1)
gene_coding_id$V2 = as.character(gene_coding_id$V2)
gene_coding_id$V2 = substr(gene_coding_id$V2,1,18)

gene_coding_id$V2 = gsub(";","",gene_coding_id$V2)
gene_coding_id_unique = data.frame(unique(gene_coding_id))
gene_coding_id_unique = data.frame(gene_coding_id_unique[,-1])
colnames(gene_coding_id_unique) = c("geneID")

id_change = read.table("rzli/ensembl_pb_gene_id",sep = '\t')
id_change$V1 = as.character(id_change$V1)
id_change$V2 = as.character(id_change$V2)

colnames(id_change) = c('en_id','geneID')

data_08 = read.table("rzli/tmp.txt",sep = '\t',header = T,row.names = 1)
data_08 = data_08[-c(34120:34124),-c(9,18,35)]
colnames(data_08) = c(paste(rep("DMSO_F"),1:4,sep = "_"),paste(rep("DMSO_M"),1:4,sep = "_"),
                      paste(rep("EME2_F"),1:4,sep = "_"),paste(rep("EME2_M"),1:4,sep = "_"),
                      paste(rep("EM_F"),1:4,sep = "_"),paste(rep("EM_M"),1:4,sep = "_"),
                      paste(rep("E2_F"),1:4,sep = "_"),paste(rep("E2_M"),1:4,sep = "_"))


data_64 = read.table("rzli/tmp_64.txt",sep = '\t',header = T,row.names = 1)
data_64 = data_64[-c(34120:34124),]
colnames(data_64) = c(paste(rep("DMSO_F"),1:4,sep = "_"),paste(rep("DMSO_M"),1:4,sep = "_"),
                      paste(rep("EME2_F"),1:4,sep = "_"),paste(rep("EME2_M"),1:4,sep = "_"),
                      paste(rep("EM_F"),1:4,sep = "_"),paste(rep("EM_M"),1:4,sep = "_"),
                      paste(rep("E2_F"),1:4,sep = "_"),paste(rep("E2_M"),1:4,sep = "_"))


data = data.frame(data_08[,17:20],data_64[,17:20])
time = factor(c(rep("8",4),rep("64",4)))
dds <- DESeqDataSetFromMatrix(data, DataFrame(time), design= ~ time)
dds <- DESeq(dds)
resdata <- as.data.frame(counts(dds, normalized=TRUE))
res <- results(dds)
res = res[order(res$padj),]
res = as.data.frame(res)

#row_data = data.frame(row.names(res))
#colnames(row_data) = c("geneID")

res_up=res[res$padj<0.05 & res$log2FoldChange>1,]
res_down=res[res$padj<0.05 & res$log2FoldChange<(-1),]
res_up=na.omit(res_up)
res_down=na.omit(res_down)

row_data_up = data.frame(row.names(res_up))
colnames(row_data_up) = c("geneID")
row_data_down = data.frame(row.names(res_down))
colnames(row_data_down) = c("geneID")

row_data_up_p = data.frame(row_data_up[which(substr(row_data_up[,1],1,1) == "P"),])
row_data_up_e = data.frame(row_data_up[which(substr(row_data_up[,1],1,1) == "E"),])
colnames(row_data_up_p) = c("geneID")
colnames(row_data_up_e) = c("geneID")
row_data_down_p = data.frame(row_data_down[which(substr(row_data_down[,1],1,1) == "P"),])
row_data_down_e = data.frame(row_data_down[which(substr(row_data_down[,1],1,1) == "E"),])
colnames(row_data_down_p) = c("geneID")
colnames(row_data_down_e) = c("geneID")

new_e_up = merge(row_data_up_p,id_change,by="geneID")
new_e_down = merge(row_data_down_p,id_change,by="geneID")

en_id_up = data.frame(new_e_up$en_id)
colnames(en_id_up) = c("geneID")
new_up = rbind(en_id_up,row_data_up_e)
en_id_down = data.frame(new_e_down$en_id)
colnames(en_id_down) = c("geneID")
new_down = rbind(en_id_down,row_data_down_e)

new_lnc_up = merge(gene_lnc_id_unique,row_data_up,by="geneID",all = F)
new_lnc_down = merge(gene_lnc_id_unique,row_data_down,by="geneID",all = F)
new_coding_up = merge(gene_coding_id_unique,row_data_up,by="geneID",all = F)
new_coding_down = merge(gene_coding_id_unique,row_data_down,by="geneID",all = F)
  
library(DOSE)
library(org.Dr.eg.db)
#library(topGO)
library(clusterProfiler)
#library(pathview)
library(AnnotationHub)

BP_up <- enrichGO(new_up$geneID,
                  "org.Dr.eg.db",
                  keyType="ENSEMBL",
                  ont="BP"
)

barplot(BP_up, showCategory=20,title = "Biological process")

BP_down <- enrichGO(new_down$geneID,
                    "org.Dr.eg.db",
                    keyType="ENSEMBL",
                    ont="BP"
)

barplot(BP_down, showCategory=20,title = "Biological process")

grep("ENSDARG00000117298",row_data_up$geneID)
grep("ENSDARG00000117298",row_data_down$geneID)

library(edgeR)
library(limma)
group = factor(c(rep("8",4),rep("64",4)))
design <- model.matrix(~0+group)
rownames(design) = colnames(data)
colnames(design) <- levels(group)
DGElist <- DGEList(counts = data, group = group)
DGElist <- calcNormFactors( DGElist )
v <- voom(DGElist, design, plot = TRUE)

dd = voom(data, design,plot = T)


#####################3
###################
######################
library(DESeq2)

data = read.table("/home/ug0007/rzli/tmp.txt",sep = '\t',header = T,row.names = 1)
data = data[-c(34120:34124),-53]
colnames(data) = c(paste(rep("Null_0_F"),1:2,sep = "_"),paste(rep("Null_0_M"),1:2,sep = "_"),
                   paste(rep("DMSO_8_F"),1:2,sep = "_"),paste(rep("DMSO_8_M"),1:2,sep = "_"),
                   paste(rep("EME2_8_F"),1:2,sep = "_"),paste(rep("EME2_8_M"),1:2,sep = "_"),
                   paste(rep("EM_8_F"),1:2,sep = "_"),paste(rep("EM_8_M"),1:2,sep = "_"),
                   paste(rep("E2_8_F"),1:2,sep = "_"),paste(rep("E2_8_M"),1:2,sep = "_"),
                   paste(rep("DMSO_64_F"),1:2,sep = "_"),paste(rep("DMSO_64_M"),1:2,sep = "_"),
                   paste(rep("EME2_64_F"),1:2,sep = "_"),paste(rep("EME2_64_M"),1:2,sep = "_"),
                   paste(rep("EM_64_F"),1:2,sep = "_"),paste(rep("EM_64_M"),1:2,sep = "_"),
                   paste(rep("E2_64_F"),1:2,sep = "_"),paste(rep("E2_64_M"),1:2,sep = "_"),
                   paste(rep("DMSO_96_F"),1:2,sep = "_"),paste(rep("DMSO_96_M"),1:2,sep = "_"),
                   paste(rep("EME2_96_F"),1:2,sep = "_"),paste(rep("EME2_96_M"),1:2,sep = "_"),
                   paste(rep("EM_96_F"),1:2,sep = "_"),paste(rep("EM_96_M"),1:2,sep = "_"),
                   paste(rep("E2_96_F"),1:2,sep = "_"),paste(rep("E2_96_M"),1:2,sep = "_")
)

data_08 = data[,c(5:20)]
data_64 = data[,c(21:36)]
data_96 = data[,c(37:52)]

cpm = apply(data,2,function(x) {x/sum(x)*1000000})

treat = factor(c(rep("D",4),rep("EME2",4),rep("EM",4),rep("E2",4)))
sex = factor(c(rep("f",2),rep("m",2),rep("f",2),rep("m",2),rep("f",2),rep("m",2),rep("f",2),rep("m",2)))

dds <- DESeqDataSetFromMatrix(data_64, DataFrame(treat,sex), design= ~ sex + treat )
dds <- DESeq(dds)

resdata <- as.data.frame(counts(dds, normalized=TRUE))
#res <- results(dds,contrast = list(c("treat_EM_vs_D")))
#res <- results(dds,contrast = list(c("sex_m_vs_f","sexm.treatEM")))
res <- results(dds,contrast = list(c("sex_m_vs_f")))
res = res[order(res$padj),]
res = as.data.frame(res)

res_up=res[res$padj<0.05 & res$log2FoldChange>1,]
res_down=res[res$padj<0.05 & res$log2FoldChange<(-1),]
res_up=na.omit(res_up)
res_down=na.omit(res_down)
res_up$geneID = row.names(res_up)
res_down$geneID = row.names(res_down)

#######
library(DOSE)
library(org.Dr.eg.db)
#library(topGO)
library(clusterProfiler)
#library(pathview)
library(AnnotationHub)

BP_up <- enrichGO(res_up$geneID,
                  "org.Dr.eg.db",
                  keyType="ENSEMBL",
                  ont="BP"
)

barplot(BP_up, showCategory=20,title = "Biological process")

BP_down <- enrichGO(res_down$geneID,
                    "org.Dr.eg.db",
                    keyType="ENSEMBL",
                    ont="BP"
)

barplot(BP_down, showCategory=20,title = "Biological process")

up = BP_up@result
down = BP_down@result
write.csv(up,"/home/ug0007/rzli/bio/64_DMSO_MvsDMSO_F_up.csv")
write.csv(down,"/home/ug0007/rzli/bio/64_DMSO_MvsDMSO_F_down.csv")
############

library(pheatmap)
a = pheatmap(cor(resdata))

#########################
#########################
aa_8 = read.table("/home/ug0007/rzli/gene_fpkm_8",header = T)
ab_8 = data.frame(aa_8$SEX0108.1.1_FRRB19H001013.1a_quant_0,
                  aa_8$SEX0108.1.2_FRRB19H001014.1a_quant_0,
                  aa_8$SEX0108.2.1_FRRB19H001029.1a_quant_0,
                  aa_8$SEX0108.2.2_FRRB19H001030.1a_quant_0,
                  aa_8$SEX0208.1.1_FRRB19H001017.1a_quant_0,
                  aa_8$SEX0208.1.2_FRRB19H001018.1a_quant_0,
                  aa_8$SEX0208.2.1_FRRB19H001033.1a_quant_0,
                  aa_8$SEX0208.2.2_FRRB19H001034.1a_quant_0,
                  aa_8$SEX0308.1.1_FRRB19H001021.1a_quant_0,
                  aa_8$SEX0308.1.2_FRRB19H001022.1a_quant_0,
                  aa_8$SEX0308.2.1_FRRB19H001037.1a_quant_0,
                  aa_8$SEX0308.2.2_FRRB19H001038.1a_quant_0,
                  aa_8$SEX0408.1.1_FRRB19H001025.1a_quant_0,
                  aa_8$SEX0408.1.2_FRRB19H001026.1a_quant_0,
                  aa_8$SEX0408.2.1_FRRB19H001041.1a_quant_0,
                  aa_8$SEX0408.2.2_FRRB19H001042.1a_quant_0)

row.names(ab_8) = aa_8$tracking_id
colnames(ab_8) = c(paste(rep("DMSO_8_F"),1:2,sep = "_"),paste(rep("DMSO_8_M"),1:2,sep = "_"),
                   paste(rep("EME2_8_F"),1:2,sep = "_"),paste(rep("EME2_8_M"),1:2,sep = "_"),
                   paste(rep("EM_8_F"),1:2,sep = "_"),paste(rep("EM_8_M"),1:2,sep = "_"),
                   paste(rep("E2_8_F"),1:2,sep = "_"),paste(rep("E2_8_M"),1:2,sep = "_"))

a = pheatmap(cor(ab_8))


aa_64 = read.table("/home/ug0007/rzli/gene_fpkm_64",header = T)
ab_64 = data.frame(aa_64$SEX0164.1.1_FRRB19H000825.1a_quant_0,
                  aa_64$SEX0164.1.2_FRRB19H001045.1a_quant_0,
                  aa_64$SEX0164.2.1_FRRB19H000975.1a_quant_0,
                  aa_64$SEX0164.2.2_FRRB19H000976.1a_quant_0,
                  aa_64$SEX0264.1.1_FRRB19H000826.1a_quant_0,
                  aa_64$SEX0264.1.2_FRRB19H001048.1a_quant_0,
                  aa_64$SEX0264.2.1_FRRB19H000979.1a_quant_0,
                  aa_64$SEX0264.2.2_FRRB19H000980.1a_quant_0,
                  aa_64$SEX0364.1.1_FRRB19H000967.1a_quant_0,
                  aa_64$SEX0364.1.2_FRRB19H000968.1a_quant_0,
                  aa_64$SEX0364.2.1_FRRB19H000983.1a_quant_0,
                  aa_64$SEX0364.2.2_FRRB19H000984.1a_quant_0,
                  aa_64$SEX0464.1.1_FRRB19H000971.1a_quant_0,
                  aa_64$SEX0464.1.2_FRRB19H000972.1a_quant_0,
                  aa_64$SEX0464.2.1_FRRB19H000987.1a_quant_0,
                  aa_64$SEX0464.2.2_FRRB19H000988.1a_quant_0)

row.names(ab_64) = aa_64$tracking_id
colnames(ab_64) = c(paste(rep("DMSO_64_F"),1:2,sep = "_"),paste(rep("DMSO_64_M"),1:2,sep = "_"),
                   paste(rep("EME2_64_F"),1:2,sep = "_"),paste(rep("EME2_64_M"),1:2,sep = "_"),
                   paste(rep("EM_64_F"),1:2,sep = "_"),paste(rep("EM_64_M"),1:2,sep = "_"),
                   paste(rep("E2_64_F"),1:2,sep = "_"),paste(rep("E2_64_M"),1:2,sep = "_"))


a = pheatmap(cor(ab_64))


aa_96 = read.table("/home/ug0007/rzli/gene_fpkm_96",header = T)
ab_96 = data.frame(aa_96$SEX0196.1.1_FRRB19H000667.1a_quant_0,
                   aa_96$SEX0196.1.2_FRRB19H000671.1a_quant_0,
                   aa_96$SEX0196.2.1_FRRB19H000675.1a_quant_0,
                   aa_96$SEX0196.2.2_FRRB19H000679.1a_quant_0,
                   aa_96$SEX0296.1.1_FRRB19H000668.1a_quant_0,
                   aa_96$SEX0296.1.2_FRRB19H000672.1a_quant_0,
                   aa_96$SEX0296.2.1_FRRB19H000676.1a_quant_0,
                   aa_96$SEX0296.2.2_FRRB19H000680.1a_quant_0,
                   aa_96$SEX0396.1.1_FRRB19H000669.1a_quant_0,
                   aa_96$SEX0396.1.2_FRRB19H000673.1a_quant_0,
                   aa_96$SEX0396.2.1_FRRB19H000677.1a_quant_0,
                   aa_96$SEX0396.2.2_FRRB19H000681.1a_quant_0,
                   aa_96$SEX0496.1.1_FRRB19H000670.1a_quant_0,
                   aa_96$SEX0496.1.2_FRRB19H000674.1a_quant_0,
                   aa_96$SEX0496.2.1_FRRB19H000678.1a_quant_0,
                   aa_96$SEX0496.2.2_FRRB19H000682.1a_quant_0)

row.names(ab_96) = aa_96$tracking_id
colnames(ab_96) = c(paste(rep("DMSO_96_F"),1:2,sep = "_"),paste(rep("DMSO_96_M"),1:2,sep = "_"),
                    paste(rep("EME2_96_F"),1:2,sep = "_"),paste(rep("EME2_96_M"),1:2,sep = "_"),
                    paste(rep("EM_96_F"),1:2,sep = "_"),paste(rep("EM_96_M"),1:2,sep = "_"),
                    paste(rep("E2_96_F"),1:2,sep = "_"),paste(rep("E2_96_M"),1:2,sep = "_"))


a = pheatmap(cor(ab_96))

ab_96['ENSDARG00000041348',]
ab_96['ENSDARG00000037491',]

ab_96['ENSDARG00000014357',]
ab_96['ENSDARG00000043923',]


#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
library(DESeq2)

data = read.table("/home/ug0007/rzli/tmp.zebrafish.txt",sep = '\t',header = T,row.names = 1)
data = data[-c(32521:32525),-73]
colnames(data) = c(
  paste(rep("DMSO_8_F"),1:4,sep = "_"),paste(rep("DMSO_8_M"),1:4,sep = "_"),
  paste(rep("EME2_8_F"),1:4,sep = "_"),paste(rep("EME2_8_M"),1:3,sep = "_"),
  paste(rep("EM_8_F"),1:3,sep = "_"),paste(rep("EM_8_M"),1:3,sep = "_"),
  paste(rep("E2_8_F"),1:3,sep = "_"),paste(rep("E2_8_M"),1:4,sep = "_"),
  paste(rep("DMSO_64_F"),1:4,sep = "_"),paste(rep("DMSO_64_M"),1:3,sep = "_"),
  paste(rep("EME2_64_F"),1:4,sep = "_"),paste(rep("EME2_64_M"),1:3,sep = "_"),
  paste(rep("EM_64_F"),1:3,sep = "_"),paste(rep("EM_64_M"),1:3,sep = "_"),
  paste(rep("E2_64_M"),1:2,sep = "_"),
  paste(rep("DMSO_96_F"),1:4,sep = "_"),paste(rep("DMSO_96_M"),1:2,sep = "_"),
  paste(rep("EME2_96_F"),1:4,sep = "_"),paste(rep("EME2_96_M"),1:2,sep = "_"),
  paste(rep("EM_96_F"),1:4,sep = "_"),paste(rep("EM_96_M"),1:2,sep = "_"),
  paste(rep("E2_96_F"),1:4,sep = "_"))

data_08 = data[,1:28]
data_64 = data[,29:50]
data_96 = data[,51:72]

##########################8
treat = factor(c(rep("DMSO",8),rep("EME2",7),rep("EM",6),rep("E2",7)))
sex = factor(c(rep("F",4),rep("M",4),rep("F",4),rep("M",3),rep("F",3),rep("M",3),rep("F",3),rep("M",4)))

dds <- DESeqDataSetFromMatrix(data_08, DataFrame(treat,sex), design= ~ sex + treat + sex:treat)
dds <- DESeq(dds)

resdata <- as.data.frame(counts(dds, normalized=TRUE))
res <- results(dds)
res = res[order(res$padj),]
res = as.data.frame(res)

res_dmso_m_f = results(dds,name = c("sex_M_vs_F"))
res_dmso_m_f = res_dmso_m_f[order(res_dmso_m_f$padj),]
res_dmso_m_f = as.data.frame(res_dmso_m_f)

res_dmso_m_f_up=res_dmso_m_f[res_dmso_m_f$padj<0.05 & res_dmso_m_f$log2FoldChange>1,]
res_dmso_m_f_down=res_dmso_m_f[res_dmso_m_f$padj<0.05 & res_dmso_m_f$log2FoldChange<(-1),]
res_dmso_m_f_up=na.omit(res_dmso_m_f_up)
res_dmso_m_f_down=na.omit(res_dmso_m_f_down)


res_em_m_f = results(dds,list(c("sex_M_vs_F","sexM.treatEM")))
res_em_m_f = res_em_m_f[order(res_em_m_f$padj),]
res_em_m_f = as.data.frame(res_em_m_f)

res_em_m_f_up=res_em_m_f[res_em_m_f$padj<0.05 & res_em_m_f$log2FoldChange>1,]
res_em_m_f_down=res_em_m_f[res_em_m_f$padj<0.05 & res_em_m_f$log2FoldChange<(-1),]
res_em_m_f_up=na.omit(res_em_m_f_up)
res_em_m_f_down=na.omit(res_em_m_f_down)

res_em_m_dmso_f = results(dds,name = c("sexM.treatEM"))
res_em_m_dmso_f = res_em_m_dmso_f[order(res_em_m_dmso_f$padj),]
res_em_m_dmso_f = as.data.frame(res_em_m_dmso_f)

res_em_m_dmso_f_up=res_em_m_dmso_f[res_em_m_dmso_f$padj<0.05 & res_em_m_dmso_f$log2FoldChange>1,]
res_em_m_dmso_f_down=res_em_m_dmso_f[res_em_m_f$padj<0.05 & res_em_m_dmso_f$log2FoldChange<(-1),]
res_em_m_dmso_f_up=na.omit(res_em_m_dmso_f_up)
res_em_m_dmso_f_down=na.omit(res_em_m_dmso_f_down)


library(pheatmap)
a = pheatmap(cor(ab_96,method = "spearman"),display_numbers = T)


############################64
treat = factor(c(rep("DMSO",7),rep("EME2",7),rep("EM",6),rep("E2",2)))
sex = factor(c(rep("F",4),rep("M",3),rep("F",4),rep("M",3),rep("F",3),rep("M",3),rep("M",2)))

dds <- DESeqDataSetFromMatrix(data_64, DataFrame(treat,sex), design= ~ sex + treat + sex:treat)
dds <- DESeq(dds)

resdata <- as.data.frame(counts(dds, normalized=TRUE))
res <- results(dds)
res = res[order(res$padj),]
res = as.data.frame(res)

res_em_dmso_m_f = results(dds,name = c("sex_M_vs_F"))

library(pheatmap)
a = pheatmap(cor(ab_96,method = "spearman"),display_numbers = T)



##############################################
library(DESeq2)
library(BiocParallel)

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


data_08 = data[,1:32]
data_64 = data[,33:64]
data_96 = data[,65:96]

data_08_use = data_08[-c(32521:34125),-c(25,17,18,16,29,32)]
data_64_use = data_64[-c(32521:34125),-c(20,7,8,29,32,13,14,21,24)]

########################
treat_08 = factor(c(rep("DMSO",8),rep("EME2",7),rep("EM",6),rep("E2",5)))
sex_08 = factor(c(rep("F",4),rep("M",4),rep("F",4),rep("M",3),rep("F",2),rep("M",4),rep("F",3),rep("M",2)))

#sex_08 <- relevel(sex_08, "M")

dds_08 <- DESeqDataSetFromMatrix(data_08_use, DataFrame(treat_08,sex_08), design= ~ sex_08 + treat_08 + sex_08:treat_08)
#dds_08 <- DESeq(dds_08,parallel = T)
dds_08 <- DESeq(dds_08)
############################################
resdata_08 <- as.data.frame(counts(dds_08, normalized=TRUE))
res_08 <- results(dds_08)
res_08 <- results(dds_08,contrast = list(c("treat_08_EM_vs_DMSO")))
#res_08 <- results(dds_08,contrast = list(c("sex_08_M_vs_F","sex_08M.treat_08EM")))
#res_08 <- results(dds_08,contrast = list(c("sex_08_M_vs_F")))

#res_08 <- results(dds_08,contrast = list(c("treat_08_EM_vs_DMSO")))

res_08 = res_08[order(res_08$padj),]
res_08 = as.data.frame(res_08)

res_08_up=res_08[res_08$padj<0.05 & res_08$log2FoldChange>1,]
res_08_down=res_08[res_08$padj<0.05 & res_08$log2FoldChange<(-1),]
res_08_up=na.omit(res_08_up)
res_08_down=na.omit(res_08_down)
res_08_up$geneID = row.names(res_08_up)
res_08_down$geneID = row.names(res_08_down)

p = pheatmap(cor(resdata_08,method = "spearman"))


library(DOSE)
library(org.Dr.eg.db)
#library(topGO)
library(clusterProfiler)
#library(pathview)
library(AnnotationHub)

BP_up <- enrichGO(res_08_up$geneID,
                  "org.Dr.eg.db",
                  keyType="ENSEMBL",
                  ont="BP"
)

barplot(BP_up, showCategory=30,title = "Biological process")

BP_down <- enrichGO(res_08_down$geneID,
                    "org.Dr.eg.db",
                    keyType="ENSEMBL",
                    ont="BP"
)

barplot(BP_down, showCategory=30,title = "Biological process")

##################################
##################################
treat_64 = factor(c(rep("DMSO",6),rep("EME2",6),rep("EM",5),rep("E2",6)))
sex_64 = factor(c(rep("F",4),rep("M",2),rep("F",4),rep("M",2),rep("F",3),rep("M",2),rep("F",4),rep("M",2)))

#sex_64 <- relevel(sex_64, "M")

dds_64 <- DESeqDataSetFromMatrix(data_64_use, DataFrame(treat_64,sex_64), design= ~ sex_64 + treat_64 + sex_64:treat_64)
#dds_64 <- DESeq(dds_64,parallel = T)
dds_64 <- DESeq(dds_64)
resdata_64 <- as.data.frame(counts(dds_64, normalized=TRUE))
#res_64 <- results(dds_64)
#res_64 <- results(dds_64,contrast = list(c("treat_64_EM_vs_DMSO")))
#res_64 <- results(dds_64,contrast = list(c("sex_64_M_vs_F","sex_64M.treat_64EM")))
res_64 <- results(dds_64,contrast = list(c("sex_64_M_vs_F")))

#res_64 <- results(dds_64,contrast = list(c("treat_64_EM_vs_DMSO")))

res_64 = res_64[order(res_64$padj),]
res_64 = as.data.frame(res_64)

res_64_up=res_64[res_64$padj<0.05 & res_64$log2FoldChange>1,]
res_64_down=res_64[res_64$padj<0.05 & res_64$log2FoldChange<(-1),]
res_64_up=na.omit(res_64_up)
res_64_down=na.omit(res_64_down)
res_64_up$geneID = row.names(res_64_up)
res_64_down$geneID = row.names(res_64_down)

vsd <- vst(dds_64, blind = FALSE)
res_64 = assay(vsd)


p = pheatmap(cor(resdata_64,method = "spearman"))


library(DOSE)
library(org.Dr.eg.db)
#library(topGO)
library(clusterProfiler)
#library(pathview)
library(AnnotationHub)

BP_up <- enrichGO(res_64_up$geneID,
                  "org.Dr.eg.db",
                  keyType="ENSEMBL",
                  ont="BP"
)

barplot(BP_up, showCategory=30,title = "Biological process")

BP_down <- enrichGO(res_64_down$geneID,
                    "org.Dr.eg.db",
                    keyType="ENSEMBL",
                    ont="BP"
)

barplot(BP_down, showCategory=30,title = "Biological process")

##################################gesa
library(clusterProfiler)
library(dplyr)

res_08_all = rbind(res_08_up,res_08_down)

res_08_new = res_08_all[order(-res_08_all$log2FoldChange),]

res_08_eg = bitr(res_08_new$geneID,fromType = "ENSEMBL",toType = "ENTREZID",
                 OrgDb = "org.Dr.eg.db")          
comb = inner_join(res_08_eg,res_08_new,by = c("ENSEMBL"='geneID'))
genelist = comb$log2FoldChange
names(genelist) = c(comb$ENTREZID)

kk = gseKEGG(geneList = genelist,
             organism = "dre",
             nPerm=10000,
             pvalueCutoff = 0.25)
ridgeplot(kk)

gseaplot(kk, geneSetID = 1, by = "runningScore", title = kk$Description[1])


##################################################

res_64_all = rbind(res_64_up,res_64_down)

res_64_new = res_64_all[order(-res_64_all$log2FoldChange),]

res_64_eg = bitr(res_64_new$geneID,fromType = "ENSEMBL",toType = "ENTREZID",
                      OrgDb = "org.Dr.eg.db")          
comb = inner_join(res_64_eg,res_64_new,by = c("ENSEMBL"='geneID'))
genelist = comb$log2FoldChange
names(genelist) = c(comb$ENTREZID)

kk = gseKEGG(geneList = genelist,
             organism = "dre",
             nPerm=10000,
             pvalueCutoff = 0.25)
ridgeplot(kk)

gseaplot(kk, geneSetID = 1, by = "runningScore", title = kk$Description[1])

######################
######################

res_64_new = res_64[order(-res_64$log2FoldChange),]

res_64_eg = bitr(res_64_new$geneID,fromType = "ENSEMBL",toType = "ENTREZID",
                 OrgDb = "org.Dr.eg.db")          
comb = inner_join(res_64_eg,res_64_new,by = c("ENSEMBL"='geneID'))
genelist = comb$log2FoldChange
names(genelist) = c(comb$ENTREZID)

kk = gseKEGG(geneList = genelist,
             organism = "dre",
             nPerm=10000,
             pvalueCutoff = 0.25)
ridgeplot(kk)

gseaplot(kk, geneSetID = 1, by = "runningScore", title = kk$Description[2])


########################################lav
library(DESeq2)

data_DMSO = read.table("/home/ug0007/rzli/fish/lav_DMSO.txt",sep = '\t',header = T,row.names = 1)
data_DMSO = data_DMSO[-c(32521:32525),-17]
data_EM = read.table("/home/ug0007/rzli/fish/lav_EM.txt",sep = '\t',header = T,row.names = 1)
data_EM = data_EM[-c(32521:32525),-18]

data_lav = cbind(data_DMSO,data_EM)

day_lav = factor(c(rep("2",2),rep("8",3),rep("12",3),rep("16",3),rep("32",3),rep("64",2),
                   rep("2",3),rep("8",3),rep("12",3),rep("16",3),rep("32",3),rep("64",2)))
treat_lav = factor(c(rep("DMSO",16),rep("EM",17)))

dds_lav <- DESeqDataSetFromMatrix(data_lav, DataFrame(treat_lav,day_lav), design= ~ day_lav + treat_lav + day_lav:treat_lav)
dds_lav <- DESeq(dds_lav)

vsd <- vst(dds_lav, blind = FALSE)
res_lav = assay(vsd)

#res_lav <- results(dds_lav)
#res_lav <- results(dds_lav,contrast = list(c("treat_lav_EM_vs_DMSO")))
#res_lav <- results(dds_lav,contrast = list(c("sex_lav_M_vs_F","sex_lavM.treat_lavEM")))
res_lav <- results(dds_lav,contrast = list(c("sex_lav_M_vs_F")))

#res_lav <- results(dds_lav,contrast = list(c("treat_lav_EM_vs_DMSO")))

res_lav = res_lav[order(res_lav$padj),]
res_lav = as.data.frame(res_lav)

res_lav_up=res_lav[res_lav$padj<0.05 & res_lav$log2FoldChange>1,]
res_lav_down=res_lav[res_lav$padj<0.05 & res_lav$log2FoldChange<(-1),]
res_lav_up=na.omit(res_lav_up)
res_lav_down=na.omit(res_lav_down)
res_lav_up$geneID = row.names(res_lav_up)
res_lav_down$geneID = row.names(res_lav_down)

pheatmap(cor(res_lav,method = "spearman"),display_numbers = T)



#######################################32
library(DESeq2)

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


vsd <- vst(dds_32, blind = FALSE)
res_32 = assay(vsd)

pheatmap(cor(res_32,method = "spearman"),display_numbers = T)

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



#library(DOSE)
library(org.Dr.eg.db)
#library(topGO)
library(clusterProfiler)
#library(pathview)
#library(AnnotationHub)

BP_up <- enrichGO(res_32_up$geneID,
                  "org.Dr.eg.db",
                  keyType="ENSEMBL",
                  ont="BP"
)

barplot(BP_up, showCategory=30,title = "Biological process")

BP_down <- enrichGO(res_32_down$geneID,
                    "org.Dr.eg.db",
                    keyType="ENSEMBL",
                    ont="BP"
)

library(ggplot2)

p = barplot(BP_down, showCategory=20,title = "Biological process")
p + scale_x_discrete(labels = function(x)str_wrap(x,width = 30) )
  



