library(ggrepel)
library(dplyr)
library(ggplot2)


data$color <- ifelse(data$padj<0.05 & abs(data$log2FoldChange)>= 1,ifelse(data$log2FoldChange > 1,'red','blue'),'gray')
color <- c(red = "red",gray = "gray",blue = "blue")
p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$padj), 
                colour=color,
                label = data$color)) +
  geom_point(alpha=0.1, size=1.5) +
  scale_color_manual(values=c("blue", "grey","red"))+
  xlim(c(-7, 15)) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.000001),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)",
       title="Differential metabolites")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())


data$label=ifelse(data$SYMBOL == "amh" | data$SYMBOL == "vtg6" | data$SYMBOL == "vtg4" ,data$SYMBOL,"")
p+geom_text_repel(data = data, aes(x = data$log2FoldChange, 
                                   y = -log10(data$padj), 
                                   label = label),
                  size = 5,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE)
