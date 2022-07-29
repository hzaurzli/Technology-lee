#PCA tSNE in single cell
library(pheatmap)
library(Rtsne)
library(ggplot2)
library(ggfortify)
library(mvtnorm)

#Establishing normal random matrix

if (pca) {
  ng = 50
  nc = 20
  a1 = rnorm(ng*nc);dim(a1) = c(ng,nc) #random matrix a1
  a2 = rnorm(ng*nc);dim(a2) = c(ng,nc) #random matrix a2
  a3 = cbind(a1,a2)
  colnames(a3) = c(paste0("cell_01_",1:nc),paste0("cell_02_",1:nc)) #cell name
  rownames(a3)=paste("gene_",1:ng,sep = "") #gene name
  pheatmap(a3) #random matrix heatmap
  #dist(a3) #calculate distance of genes
  a3=t(a3) #calculate distance of cells
  dim(a3)
  #pca_dot = prcomp(a3,scale. = T) #pca deposition
  #p = autoplot(pca_dot)+theme_classic()+ggtitle("PCA")
  #p1=screeplot(pca_dot) #screeplot 
  df = cbind(as.data.frame(a3),group=c(rep('b1',20),rep('b2',20))) #b1,b2 distinguish color
  autoplot(prcomp(df[,1:(ncol(df)-1)]),data = df,colour = 'group')+theme_bw()
}

if(tsne){
  set.seed(42)
  tsne_out = Rtsne(a3,pca = F,perplexity = 10,theta = 0.0) #run tsne
  tsnes = tsne_out$Y
  colnames(tsnes)=c("tsne1","tsne2")
  group=c(rep('b1',20),rep('b2',20))
  ggplot(as.data.frame(tsnes),aes(x = tsne1,y = tsne2))+geom_point(aes(col=group))#tsne plot
}

