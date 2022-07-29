library(AnnotationHub)
library(clusterProfiler)
library(org.Osativa.eg.db)
org <-org.Osativa.eg.db
rice_kegg <- clusterProfiler::download_KEGG("dosa")

zmj02 = read.delim("/home/lirz/rice/Galaxy102-[DEG_analysis_on_DESeq2_gene_result].tabular",sep = '\t',header = T)
zmj02 = zmj02[zmj02$P.adj<0.05,]

ego_BP <- enrichGO(zmj02$GeneID,
                   OrgDb = org,
                   keyType = "GID",
                   pAdjustMethod = "BH",
                   ont="BP")

p1 <- dotplot(ego_BP,title = "Biological process")

ego_MF <- enrichGO(zmj02$GeneID,
                   OrgDb = org,
                   keyType = "GID",
                   pAdjustMethod = "BH",
                   ont="MF")

p2 <- dotplot(ego_MF,title = "Molecular Function")

ego_CC <- enrichGO(zmj02$GeneID,
                      OrgDb = org,
                      keyType = "GID",
                      pAdjustMethod = "BH",
                      ont="CC")

p3 <- dotplot(ego_CC,title = "Cellular Component")

#######select
BP = ego_BP@result
MF = ego_MF@result
CC = ego_CC@result

GO_gene_BP <- strsplit(BP$geneID,split="/")
GO_gene_MF <- strsplit(MF$geneID,split="/")
GO_gene_CC <- strsplit(CC$geneID,split="/")

select_GO = function(x){
  d = vector()
  for (i in 1:length(x)) {
    a = as.data.frame(x[i])
    b = as.data.frame(rep(BP[i,2],length(x[[i]])))
    c = cbind(a,b)
    colnames(c) = c("GeneID","function")
    d = rbind(d,c)
  }
  return(d)
}
 
BP_function = select_GO(GO_gene_BP)
MF_function = select_GO(GO_gene_MF)
CC_function = select_GO(GO_gene_CC)

#########up
zmj02_up=zmj02[zmj02$P.adj<0.05 & zmj02$log2.FC.> 1,]##change
zmj02_up<-na.omit(zmj02_up)

rap_id_up=as.vector(zmj02_up$GeneID)
map_id_up <- AnnotationDbi::select(org, keys = rap_id_up, 
                                   columns=c("RAP"), keytype = "GID")
map_id_up<-na.omit(map_id_up)


rap_id_up=paste0(map_id_up$RAP[!is.na(map_id_up$RAP)],"-01")
rap_id_up <- gsub("g","t",rap_id_up)

##kegg enrichment
kk_up <- enrichKEGG(rap_id_up, organism="dosa",
                    keyType = "kegg",                 
                    pvalueCutoff=0.05, pAdjustMethod="BH",                  
                    qvalueCutoff=0.1)

kk = dotplot(kk_up,title = "KEGG pathway up_regulation")

write.csv(kk_up@result,"/home/lirz/rice/table/kegg_up.csv")

########down
zmj02_down=zmj02[zmj02$P.adj<0.05 & zmj02$log2.FC.< c(-1),]##change
zmj02_down<-na.omit(zmj02_down)

rap_id_down=as.vector(zmj02_down$GeneID)
map_id_down <- AnnotationDbi::select(org, keys = rap_id_down, 
                                     columns=c("RAP"), keytype = "GID")
map_id_down<-na.omit(map_id_down)


rap_id_down=paste0(map_id_down$RAP[!is.na(map_id_down$RAP)],"-01")
rap_id_down <- gsub("g","t",rap_id_down)


##kegg enrichment
gg_down <- enrichKEGG(rap_id_down, organism="dosa",
                      keyType = "kegg",                 
                      pvalueCutoff=0.05, pAdjustMethod="BH",                  
                      qvalueCutoff=0.1)

gg = dotplot(gg_down,title = "KEGG pathway down_regulation")

write.csv(gg_down@result,"/home/lirz/rice/table/kegg_down.csv")

############summarise  
#######up
zmj02_up=zmj02[zmj02$P.adj<0.05 & zmj02$log2.FC.> 1,]##change
zmj02_up<-na.omit(zmj02_up)

zmj02_up_BP = merge(zmj02_up,BP_function,by = "GeneID")
zmj02_up_MF = merge(zmj02_up,MF_function,by = "GeneID")
zmj02_up_CC = merge(zmj02_up,CC_function,by = "GeneID")

write.csv(zmj02_up,"/home/lirz/rice/table/zmj02_up.csv")
write.csv(zmj02_up_BP,"/home/lirz/rice/table/zmj02_up_BP.csv")
write.csv(zmj02_up_MF,"/home/lirz/rice/table/zmj02_up_MF.csv")
write.csv(zmj02_up_CC,"/home/lirz/rice/table/zmj02_up_CC.csv")

######down
zmj02_down=zmj02[zmj02$P.adj<0.05 & zmj02$log2.FC.< c(-1),]##change
zmj02_down<-na.omit(zmj02_down)

zmj02_down_BP = merge(zmj02_down,BP_function,by = "GeneID")
zmj02_down_MF = merge(zmj02_down,MF_function,by = "GeneID")
zmj02_down_CC = merge(zmj02_down,CC_function,by = "GeneID")

write.csv(zmj02_down,"/home/lirz/rice/table/zmj02_down.csv")
write.csv(zmj02_down_BP,"/home/lirz/rice/table/zmj02_down_BP.csv")
write.csv(zmj02_down_MF,"/home/lirz/rice/table/zmj02_down_MF.csv")
write.csv(zmj02_down_CC,"/home/lirz/rice/table/zmj02_down_CC.csv")

###########heatmap
zmj02_TPM = read.delim("/home/lirz/rice/Galaxy103-[DEG_analysis_on_gene_TPM].tabular",sep = '\t',header = T)
row.names(zmj02_TPM) = zmj02_TPM[,1]
zmj02_TPM = zmj02_TPM[,-c(1,5)]
zmj02_TPM = zmj02_TPM[rowSums(zmj02_TPM)!=0,]


zmj02_order = rbind(zmj02_up,zmj02_down)
zmj02_order_genelist = as.data.frame(zmj02_order$GeneID)
colnames(zmj02_order_genelist) = c("GeneID")

GeneID = row.names(zmj02_TPM)
zmj02_TPM = cbind(GeneID,zmj02_TPM)

zmj02_heatmap = merge(zmj02_order_genelist,zmj02_TPM,by = 'GeneID')
row.names(zmj02_heatmap) = zmj02_heatmap[,1]
zmj02_heatmap = zmj02_heatmap[,-1]
colnames(zmj02_heatmap) = c("wt_a","wt_b","wt_c","x_1a","x_1b","x_1c")


pheatmap::pheatmap(cor(zmj02_heatmap,method = "spearman"),display_numbers = TRUE,number_color = "black")



######################################################################
######################################################################
library(AnnotationHub)
library(clusterProfiler)
library(org.Osativa.eg.db)
org <-org.Osativa.eg.db
rice_kegg <- clusterProfiler::download_KEGG("dosa")

zmj02 = read.delim("/home/lirz/rice/Galaxy102-[DEG_analysis_on_DESeq2_gene_result].tabular",sep = '\t',header = T)
zmj02 = zmj02[zmj02$P.adj<0.05,]

ego_BP <- enrichGO(zmj02$GeneID,
                   OrgDb = org,
                   keyType = "GID",
                   pAdjustMethod = "BH",
                   ont="BP")

p1 <- dotplot(ego_BP,title = "Biological process")

ego_MF <- enrichGO(zmj02$GeneID,
                   OrgDb = org,
                   keyType = "GID",
                   pAdjustMethod = "BH",
                   ont="MF")

p2 <- dotplot(ego_MF,title = "Molecular Function")

ego_CC <- enrichGO(zmj02$GeneID,
                   OrgDb = org,
                   keyType = "GID",
                   pAdjustMethod = "BH",
                   ont="CC")

p3 <- dotplot(ego_CC,title = "Cellular Component")

#######select
BP = ego_BP@result
MF = ego_MF@result
CC = ego_CC@result

GO_gene_BP <- strsplit(BP$geneID,split="/")
GO_gene_MF <- strsplit(MF$geneID,split="/")
GO_gene_CC <- strsplit(CC$geneID,split="/")

select_GO = function(x){
  d = vector()
  for (i in 1:length(x)) {
    a = as.data.frame(x[i])
    b = as.data.frame(rep(BP[i,2],length(x[[i]])))
    c = cbind(a,b)
    colnames(c) = c("GeneID","function")
    d = rbind(d,c)
  }
  return(d)
}

BP_function = select_GO(GO_gene_BP)
MF_function = select_GO(GO_gene_MF)
CC_function = select_GO(GO_gene_CC)

#########up
zmj02_up=zmj02[zmj02$P.adj<0.05 & zmj02$log2.FC.> 1,]##change
zmj02_up<-na.omit(zmj02_up)

rap_id_up=as.vector(zmj02_up$GeneID)
map_id_up <- AnnotationDbi::select(org, keys = rap_id_up, 
                                   columns=c("RAP"), keytype = "GID")
map_id_up<-na.omit(map_id_up)


rap_id_up=paste0(map_id_up$RAP[!is.na(map_id_up$RAP)],"-01")
rap_id_up <- gsub("g","t",rap_id_up)

##kegg enrichment
kk_up <- enrichKEGG(rap_id_up, organism="dosa",
                    keyType = "kegg",                 
                    pvalueCutoff=0.05, pAdjustMethod="BH",                  
                    qvalueCutoff=0.1)

kk = dotplot(kk_up,title = "KEGG pathway up_regulation")

write.csv(kk_up@result,"/home/lirz/rice/table/kegg_up.csv")

########down
zmj02_down=zmj02[zmj02$P.adj<0.05 & zmj02$log2.FC.< c(-1),]##change
zmj02_down<-na.omit(zmj02_down)

rap_id_down=as.vector(zmj02_down$GeneID)
map_id_down <- AnnotationDbi::select(org, keys = rap_id_down, 
                                     columns=c("RAP"), keytype = "GID")
map_id_down<-na.omit(map_id_down)


rap_id_down=paste0(map_id_down$RAP[!is.na(map_id_down$RAP)],"-01")
rap_id_down <- gsub("g","t",rap_id_down)


##kegg enrichment
gg_down <- enrichKEGG(rap_id_down, organism="dosa",
                      keyType = "kegg",                 
                      pvalueCutoff=0.05, pAdjustMethod="BH",                  
                      qvalueCutoff=0.1)

gg = dotplot(gg_down,title = "KEGG pathway down_regulation")

write.csv(gg_down@result,"/home/lirz/rice/table/kegg_down.csv")

############summarise  
#######up
zmj02_up=zmj02[zmj02$P.adj<0.05 & zmj02$log2.FC.> 1,]##change
zmj02_up<-na.omit(zmj02_up)

zmj02_up_BP = merge(zmj02_up,BP_function,by = "GeneID")
zmj02_up_MF = merge(zmj02_up,MF_function,by = "GeneID")
zmj02_up_CC = merge(zmj02_up,CC_function,by = "GeneID")

write.csv(zmj02_up,"/home/lirz/rice/table/zmj02_up.csv")
write.csv(zmj02_up_BP,"/home/lirz/rice/table/zmj02_up_BP.csv")
write.csv(zmj02_up_MF,"/home/lirz/rice/table/zmj02_up_MF.csv")
write.csv(zmj02_up_CC,"/home/lirz/rice/table/zmj02_up_CC.csv")

######down
zmj02_down=zmj02[zmj02$P.adj<0.05 & zmj02$log2.FC.< c(-1),]##change
zmj02_down<-na.omit(zmj02_down)

zmj02_down_BP = merge(zmj02_down,BP_function,by = "GeneID")
zmj02_down_MF = merge(zmj02_down,MF_function,by = "GeneID")
zmj02_down_CC = merge(zmj02_down,CC_function,by = "GeneID")

write.csv(zmj02_down,"/home/lirz/rice/table/zmj02_down.csv")
write.csv(zmj02_down_BP,"/home/lirz/rice/table/zmj02_down_BP.csv")
write.csv(zmj02_down_MF,"/home/lirz/rice/table/zmj02_down_MF.csv")
write.csv(zmj02_down_CC,"/home/lirz/rice/table/zmj02_down_CC.csv")

###########heatmap
zmj02_TPM = read.delim("/home/lirz/rice/Galaxy103-[DEG_analysis_on_gene_TPM].tabular",sep = '\t',header = T)
row.names(zmj02_TPM) = zmj02_TPM[,1]
zmj02_TPM = zmj02_TPM[,-c(1,5)]
zmj02_TPM = zmj02_TPM[rowSums(zmj02_TPM)!=0,]


zmj02_order = rbind(zmj02_up,zmj02_down)
zmj02_order_genelist = as.data.frame(zmj02_order$GeneID)
colnames(zmj02_order_genelist) = c("GeneID")

GeneID = row.names(zmj02_TPM)
zmj02_TPM = cbind(GeneID,zmj02_TPM)

zmj02_heatmap = merge(zmj02_order_genelist,zmj02_TPM,by = 'GeneID')
row.names(zmj02_heatmap) = zmj02_heatmap[,1]
zmj02_heatmap = zmj02_heatmap[,-1]
colnames(zmj02_heatmap) = c("wt_a","wt_b","wt_c","x_1a","x_1b","x_1c")


pheatmap::pheatmap(cor(zmj02_heatmap,method = "spearman"),display_numbers = TRUE,number_color = "black")

