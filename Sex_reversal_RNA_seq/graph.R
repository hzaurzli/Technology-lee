library(tidyverse)
#Reproducible data set
test_mtcars <- mtcars %>% group_by(cyl,am, gear) %>% summarise(mean = mean(mpg))

ggplot(test_mtcars, aes(as.factor(cyl), mean, fill=as.factor(am))) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_grid(~gear) + geom_text(aes(label = round(mean, 2)), position = position_dodge(width = 0.9), vjust = -1)


library(reshape2)
head(tips)


library(ggplot2)
sp <- ggplot(tips, aes(x=total_bill, y=tip/total_bill)) + geom_point(shape=1)

# 垂直方向进行分割
sp + facet_grid(sex ~ .)

# 水平方向分割
sp + facet_grid(. ~ sex)

# 对 sex 进行垂直分割, 对 day 进行水平分割
sp + facet_grid(sex ~ day)


# 水平分割, 分为两列
sp + facet_wrap( ~ day, ncol=2)


sp + facet_grid(sex ~ day) +
  theme(strip.text.x = element_text(size=8, angle=75),
        strip.text.y = element_text(size=12, face="bold"),
        strip.background = element_rect(colour="red", fill="#CCCCFF"))


labels <- c(Female = "Women", Male = "Men")
sp + facet_grid(. ~ sex, labeller=labeller(sex = labels))
#######################
height <- c(6,5.92,5.58,5.83,32,44,2432,4,1423,65647,5,65,6564,67,3)
wei <- c(20,15,7,12,21,432,43,543,654,23,5,67,43,23,432)
marks <- c(10,20,15,11,23,34,435,234412,432,546,342,32,546,76,54)
corr = data.frame(height,wei,marks)

table=cor(corr)

library(corrplot)
corlor1 = colorRampPalette(c("red","yellow","blue"))
corrplot(table,type = "full",col = corlor1(20))

a = corlor1(20)
t = c(a,a,a,a)
corrplot(table,col = t)
##############
library_id = metadata_seq$library_id[1:5]
data_size = metadata_seq$species_datasize_G[1:5]
library_id = metadata_seq$library_id[104:108]
data_size = metadata_seq$species_datasize_G[104:108]

data1 = data.frame(library_id,data_size)
data2 = data.frame(library_id,data_size)
data = rbind(data1,data2)

data$library_id = as.character(data$library_id)
data$data_size = as.numeric(data$data_size)
for (i in 1:5) {
  data[i,1] = paste(data[i,1],"z",sep="")
}

for (i in 6:10) {
  data[i,1] = paste(data[i,1],"m",sep="")
}


library(ggplot2)
ggplot(data,aes(x = data$library_id,y = data$data_size)) + 
  geom_bar(stat="identity") +theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))



