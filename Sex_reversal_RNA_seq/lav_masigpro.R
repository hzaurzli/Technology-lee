data_D = read.table("E:/R file/lav_DMSO.txt",sep = '\t',header = T,row.names = 1)
data_D = data_D[-c(32521:32525),-17]
data_D = data_D[,-c(7,11,13)]

data_E = read.table("E:/R file/lav_EM.txt",sep = '\t',header = T,row.names = 1)
data_E = data_E[-c(32521:32525),-18]
data = cbind(data_D,data_E)

cpm = apply(data,2,function(x) {x/sum(x)*1000000})
cpm_D = apply(data_D,2,function(x) {x/sum(x)*1000000})

condition = data.frame( time = c(rep("1",2),rep("2",3),rep("3",2),rep("4",2),rep("5",2),rep("6",2),
                                 rep("1",3),rep("2",3),rep("3",3),rep("4",3),rep("5",3),rep("6",2)),
                        
                        rep = c(rep("1",2),rep("2",3),rep("3",2),rep("4",2),rep("5",2),rep("6",2),
                               rep("7",3),rep("8",3),rep("9",3),rep("10",3),rep("11",3),rep("12",2)),
                        control = c(rep("1",13),rep("0",17)),
                        EM = c(rep("0",13),rep("1",17)),
                        DMSO = c(rep("1",13),rep("0",17))
                        )


condition_D = data.frame( 
                        time = c(rep(2,2),rep(8,3),rep(12,2),rep(16,2),rep(32,2),rep(64,2)),
                        rep = c(rep("1",2),rep("2",3),rep("3",2),rep("4",2),rep("5",2),rep("6",2)),
                        #control = c(rep("1",13)),
                        DMSO = c(rep("1",13))
)


condition$time = as.numeric(condition$time)
condition$rep = as.numeric(condition$rep)
condition$DMSO = as.numeric(condition$DMSO)
condition$control = as.numeric(condition$control)
condition$DMSO = as.numeric(condition$DMSO)
condition$EM = as.numeric(condition$EM)

row.names(condition) =colnames(cpm)
condition = as.matrix(condition)

####################
condition_D$time = as.numeric(condition_D$time)
condition_D$rep = as.numeric(condition_D$rep)
#condition_D$control = as.numeric(condition_D$control)
condition_D$DMSO = as.numeric(condition_D$DMSO)

row.names(condition_D) =colnames(cpm_D)
condition_D = as.matrix(condition_D)
                 
library(maSigPro)

design <- make.design.matrix(condition_D,degree = 6)

fit <- p.vector(cpm_D, design)#, Q = 0.05, MT.adjust = "BH", min.obs = 20)
tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
sigs <- get.siggenes(tstep, rsq = 0.15, vars = "groups")

#suma2Venn(sigs$summary[, c(2:4)])
#suma2Venn(sigs$summary[, c(1:4)])

#sigs$sig.genes$DMSO
a = see.genes(sigs$sig.genes$DMSO, 
              show.fit = T, dis =design$dis, 
              cluster.method="hclust" ,
              cluster.data = 1, k = 4)


names(sigs$sig.genes$DMSO)

group = data.frame(a$cut)

cluster1 <- cpm_D[which(group=="2"),]

PlotGroups(cluster1, edesign = condition_D)

#########################
lnc_id = read.csv("E:/R file/lnc_id.csv")
lnc_id = lnc_id[,-1]
coding_id = read.csv("E:/R file/coding_id.csv")
coding_id = coding_id[,-1]

gene = data.frame(a$cut)
gene_id = row.names(gene)
gene = cbind(gene,gene_id)

ddg_lnc = merge(gene,lnc_id,by='gene_id',all = F)
ddg_coding = merge(gene,coding_id,by='gene_id',all = F)

##############
coding_lnc = data.frame(
  type = c(rep("other",158),rep("coding",10095),rep("lnc",470),rep("other",4067),rep("coding",15337),rep("lnc",2393)),
  group = c(rep("DDG",158),rep("DDG",10095),rep("DDG",470),rep("non-DDG",4067),rep("non-DDG",15337),rep("non-DDG",2393))
)

coding_lnc = data.frame(
  Gene_type = c("LncRNA","LncRNA",
                "Protein coding","Protein coding"),
  Type = c("sex_bias","sex_unbias","sex_bias","sex_unbias"),
  Gene_percent = c(1287,1576,14791,10641),
  Gene_number = c(1287,1576,14791,10641)
)

coding_lnc = data.frame(
  Gene_type = c(rep('LncRNA',470),rep('Protein coding',2393),
                rep('LncRNA',10095),rep('Protein coding',15337)),
  Type = c(rep('DDG',470),rep('non-DDG',2393),
           rep('DDG',10095),rep('non-DDG',15337))
)


coding_lnc$type = factor(coding_lnc$Gene_type, levels=c('other','lnc','coding'))

library(ggplot2)

p = ggplot(coding_lnc, aes(x = Gene_type, y = Gene_percent, fill = Type)) + 
  geom_bar(stat = "identity",width = 0.4,position = "fill") + 
  geom_text(mapping = aes(label = Gene_number),color = 'white',size = 5,position = position_fill(0.5)) +
  scale_fill_manual(values=c('#FF7F00','#386CB0')) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(
    axis.ticks.length = unit(0.5,'cm'),
    panel.grid = element_blank(),
    axis.title.x=element_text(size=15),
    axis.title.y=element_text(size=15),
    axis.text.x=element_text(size=14),
    axis.text.y=element_text(size=14),
    legend.title = element_text(size=13),
    legend.text = element_text(size = 12)
  ) +  coord_flip() 
  
p


library(ggplot2) #加载ggplot2包
library(dplyr) #加载dplyr包
library(ggstatsplot)


coding_lnc %>% 
  ggplot(aes(x = Gene_type, fill = Type)) + 
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set2') + #设置颜色板
  theme_bw() + #设置主题
  theme(panel.grid = element_blank()) +  #清除网格线
  labs(y = 'Percent') + #设置y轴名为‘Percent’
  scale_y_continuous(labels = scales::percent) + #设置y轴的标签形式为百分比
  coord_flip() #旋转坐标轴


require("RColorBrewer")
display.brewer.all(type = "qual")
paired=brewer.pal(n = 8, name = "Set1")




data<-data.frame(count=c(39,36,19,6), category=c("a","b","c","d"))
data$fraction = data$count / sum(data$count)
data = data[order(data$fraction), ]
data$ymax = cumsum(data$fraction)
data$ymin = c(0, head(data$ymax, n=-1))


fill <- c("blue3","cyan3","darkgrey","forestgreen")

library(ggplot2)

p1 = ggplot(data, aes(fill=category, ymax=ymax, ymin=ymin, xmax=4, xmin=3.5)) +
  geom_rect(colour="White") +
  coord_polar(theta="y") +
  scale_fill_manual(values=fill)+
  theme_bw()+
  geom_label(aes(label=paste(data$fraction*100,"%"),x=4,y=
                   (ymin+ymax)/2),inherit.aes = F)+
  theme(panel.grid=element_blank())+
  theme(axis.ticks=element_blank()) +     
  xlim(c(0, 4)) +
  theme(axis.text=element_blank()) +
  theme(legend.text=element_text(color=fill,size=12))+
  theme(legend.key.size=unit(2,'lines'))+
  theme(legend.key=element_rect(size=5))+
  labs(title="donut plot")




sex_bias = data.frame(
  Gene_type = c("LncRNA","LncRNA",
                "Protein coding","Protein coding"),
  Type = c("sex_bias","sex_unbias","sex_bias","sex_unbias"),
  Gene_percent = c(1303,1560,16775,8657),
  Gene_number = c(1303,1560,16775,8657)
)

library(ggplot2)

p = ggplot(sex_bias, aes(x = Gene_type, y = Gene_percent, fill = Type)) + 
  geom_bar(stat = "identity",width = 0.4,position = "fill") + 
  geom_text(mapping = aes(label = Gene_number),color = 'white',size = 5,position = position_fill(0.5)) +
  scale_fill_manual(values=c('#FF7F00','#386CB0')) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(
    axis.ticks.length = unit(0.5,'cm'),
    panel.grid = element_blank(),
    axis.title.x=element_text(size=15),
    axis.title.y=element_text(size=15),
    axis.text.x=element_text(size=14),
    axis.text.y=element_text(size=14),
    legend.title = element_text(size=13),
    legend.text = element_text(size = 12)
  ) +  coord_flip() 



