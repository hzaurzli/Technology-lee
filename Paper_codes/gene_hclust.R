# 基因 hclust 聚类并绘制聚类热图和表达趋势图

# 加载R包
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggnewscale)
library(circlize)

# 读取数据
exps <- read.delim('c:/Users/admin/Desktop/protein_exp.txt',row.names = 1,header = T)

# 查看数据
head(exps,3)

# 归一化方便热图展示
hclust_matrix <- exps %>% t() %>% scale() %>% t()

# 计算聚类并聚类
gene_hclust <- hclust(dist(hclust_matrix),method = 'complete')

# 绘制聚类树
plot(gene_hclust, labels = FALSE)
abline(h = 3.8, col = "red", lwd = 2)

# 分成12类
gene_cluster <- cutree(gene_hclust, k = 12) %>% 
  # 转换为tibble格式
  enframe() %>% 
  # 添加列名
  rename(gene = name, cluster = value)

# 查看分类
head(gene_cluster)

# 可视化cluster
my_cluster <- hclust_matrix %>% as.data.frame()
my_cluster[['gene']] <- rownames(my_cluster)

# 加入cluster信息
my_cluster <- my_cluster %>% 
  inner_join(gene_cluster, by = "gene")

# 查看内容
head(my_cluster,3)

# 查看每类基因数量
table(gene_cluster$cluster)

# wide to long
df <- melt(my_cluster,id.vars = c('gene','cluster'))

# 添加cluster name
df$cluster_name <- paste("Cluster ",df$cluster,sep = '')

# cluster 因子化
df$cluster_name <- factor(df$cluster_name,levels = paste("Cluster ",1:12,sep = ''))


# median expression in each cluster
line <- ggplot(df,aes(x = variable,y = value,group = gene)) +
  geom_line(alpha = 0.2,color = 'grey20') +
  geom_line(stat = "summary", fun = "median", colour = "brown", size = 3.5, 
            aes(group = 1)) +
  theme_bw(base_size = 16) +
  scale_x_discrete(labels = c("zygote","cell 2","cell 4","cell 8","morula","blasto")) +
  theme(axis.text.x = element_text(color = 'black'),
        axis.ticks.length = unit(0.25,'cm')) +
  xlab('') + ylab('Relative expression') +
  facet_wrap(~cluster_name,ncol = 3,scales = 'free')

line

# 热图
dp <- df %>% arrange(desc(cluster))

# gene 因子化
dp$gene <- factor(dp$gene,levels = unique(dp$gene))

# 绘制热图
heatmap <- ggplot(dp,aes(x = variable,y = gene,group = cluster_name)) +
  geom_tile(aes(fill = value),width = 1) +
  theme_bw(base_size = 16) +
  xlab('') + ylab('') +
  scale_x_discrete(position = 'top',
                   labels = c("zygote","cell 2","cell 4","cell 8","morula","blasto")) +
  scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0,
                       name = "Z-score") +
  # 添加cluster注释
  new_scale('fill') +
  geom_tile(aes(x = 0.4,y = gene,fill = cluster_name),width = 0.25) +
  scale_fill_manual(values = rand_color(12)) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 0,color = 'black',size = 16),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'left')

heatmap

library(patchwork)
# 拼图
heatmap + line + plot_layout(widths = c(0.5,1))
