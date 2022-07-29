#!/usr/bin/R
args = commandArgs(T)

count_table_path = args[1]


#install.packages("https://github.com/xuzhougeng/org.Osativa.eg.db/releases/download/v0.01/org.Osativa.eg.db.tar.gz",
                 #repos = NULL,
                 #type="source")

library(AnnotationHub)
library(clusterProfiler)
library(org.Osativa.eg.db)
org <-org.Osativa.eg.db
rice_kegg <- clusterProfiler::download_KEGG("dosa")

##Read file
preprocess1 = function(count_table_path){
  zmj103=read.table(count_table_path,sep = "\t",header = T)
  ##up
  zmj103_up=zmj103[zmj103$P.adj<0.05 & zmj103$log2.FC.> 1,]
  zmj103_up<-na.omit(zmj103_up)
  return(zmj103_up)
}


preprocess2 = function(count_table_path){
  zmj103=read.table(count_table_path,sep = "\t",header = T)
  ##down
  zmj103_down=zmj103[zmj103$P.adj<0.05 & zmj103$log2.FC.< c(-1),]
  zmj103_down<-na.omit(zmj103_down)
  return(zmj103_down)
}

###kegg
kegg1 = function(count_table_path){
  rap_id_up = preprocess1(count_table_path)$GeneID
  map_id_up <- AnnotationDbi::select(org, keys = rap_id_up,
                                columns=c("RAP"), keytype = "GID")
  map_id_up <- na.omit(map_id_up)
  rap_id_up=paste0(map_id_up$RAP[!is.na(map_id_up$RAP)],"-01")
  rap_id_up <- gsub("g","t",rap_id_up)
  ##kegg enrichment
  kk <- enrichKEGG(rap_id_up, organism="dosa",
                   keyType = "kegg",
                   pvalueCutoff=0.05, pAdjustMethod="BH",
                   qvalueCutoff=0.1)

  return(kk)
}


kegg2 = function(count_table_path){
  rap_id_down = preprocess2(count_table_path)$GeneID
  map_id_down <- AnnotationDbi::select(org, keys = rap_id_down,
                                columns=c("RAP"), keytype = "GID")
  map_id_down <- na.omit(map_id_down)
  rap_id_down=paste0(map_id_down$RAP[!is.na(map_id_down$RAP)],"-01")
  rap_id_down <- gsub("g","t",rap_id_down)
  ##kegg enrichment
  gg <- enrichKEGG(rap_id_down, organism="dosa",
                   keyType = "kegg",
                   pvalueCutoff=0.05, pAdjustMethod="BH",
                   qvalueCutoff=0.1)

  return(gg)
}


write.csv(preprocess1,"/home/lirz/zmj/data/up.csv")
write.csv(preprocess2,"/home/lirz/zmj/data/down.csv")
write.csv(kegg1,"/home/lirz/zmj/data/kegg_up.csv")
write.csv(kegg2,"/home/lirz/zmj/data/kegg_down.csv")
