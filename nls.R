ds = read.table("/home/lirz/cutoff/genes.fpkm_table",sep = '\t',header = T)
row.names(ds) = ds$tracking_id
ds =ds[,-1]

ee = c(0,0.16,0.32,0.48,0.64,0.8,1.6,2.4,3.2,4,4.8,5.6,6.4,7.2,8)
uu = select_fpkm(ds)
uu = c(0,uu)
ii = round(uu/18300,4)
data = as.data.frame(cbind(ee,uu,ii))


select_fpkm <- function(x) {
  aa = c()
  for (i in 1:14) {
    gt = x[which(x[,i] > 1),]
    num = nrow(gt)
    aa[i] = num
    
  }
  return(aa)
}


sel_fpkm <- function(x) {
  aa = c()
  for (i in 1:14) {
    gt = x[which(x[,i] > 0),]
    num = nrow(gt)
    aa[i] = num
    
  }
  return(aa)
}


f <- nls(uu ~ -a/(ee+c) + b,start=list(a = 400, b = -49, c = 0.025), data=data, trace=T)

library(ggplot2)
library(scales)
ggplot(data,aes(ee, uu))+
  geom_point(size=3)+geom_line(aes(ee,fitted(f)),col='red')+
  xlab("Data/G")+
  ylab("Gene number")+
  geom_text(aes(label = uu, vjust = 1.1 , hjust = -0.5, angle = -45),size = 3.5, show.legend = FALSE)


fp <- nls(ii ~ -a/(ee+c) + b,start=list(a = 400, b = -49, c = 0.025), data=data, trace=T)

library(scales)
ggplot(data,aes(ee, ii))+
  geom_point(size=3)+geom_line(aes(ee,fitted(fp)),col='red')+
  xlab("Data/G")+
  ylab("Gene number percent")+
  geom_text(aes(label = uu, vjust = 1.1 , hjust = -0.5,angle = -45),size = 3.5, show.legend = FALSE)+
  scale_y_continuous(labels = percent)
