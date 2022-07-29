library("doParallel")      
library("foreach")         
cl<- makeCluster(detectCores()-1)      
registerDoParallel(cl)       
getDoParWorkers()  

code=function(){
  id_change = read.table("E:/R file/adult_sex_reversal/ensembl_pb_gene_id",sep = '\t')
  id_change$V1 = as.character(id_change$V1)
  id_change$V2 = as.character(id_change$V2)
  
  res = read.csv("E:/R file/adult_sex_reversal/res.csv",header = T,row.names = 1)
  res$geneID = as.character(res$geneID)
  
  res_pb = res[which(substr(res[,7],1,1) == "P"),]
  for (i in x) {
    for(j in y){
      if(id_change[j,2] == res_pb[i,7]){
        res_pb[i,7] = id_change[j,1]
      }
    }
  }  
}    

uu <- foreach(x=1:19449, y=1:18303,.combine=rbind) %dopar% code()
stopCluster(cl) 