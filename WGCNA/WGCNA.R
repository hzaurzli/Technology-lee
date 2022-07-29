##Download the software
source("https://bioconductor.org/biocLite.R")
biocLite(c("AnnotationDbi", "impute","GO.db", "preprocessCore"))
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
install.packages(c("WGCNA", "stringr", "reshape2"), repos=site)


library(WGCNA)
library(reshape2)
library(stringr)

##Read files
fpkm=read.csv("/home/lirz/pra2/fpkm.csv",header = T)
datTraits=read.csv("/home/lirz/pra2/datTraits.csv",header = T)
row.names(fpkm)=fpkm[,1]
fpkm=fpkm[,-1]
row.names(datTraits)=datTraits[,1]
datTraits=datTraits[,-1]


RNAseq_voom <- fpkm
WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:5000],])
datExpr0 <- WGCNA_matrix 
datExpr <- datExpr0 

##Download data
library(GEOquery)
a=getGEO('GSE48213')
metadata=pData(a[[1]])[,c(2,10,12)]
datTraits = data.frame(gsm=metadata[,1],
                       cellline=trimws(sapply(as.character(metadata$characteristics_ch1),function(x) strsplit(x,":")[[1]][2])),
                       subtype=trimws(sapply(as.character(metadata$characteristics_ch1.2),function(x) strsplit(x,":")[[1]][2]))
)


##Set soft threshold,Determine the best BETA value
#Set the network construction parameter selection range and calculate the scale-free topology matrix
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
par(mfrow = c(1,2))
cex1 = 0.9
##Two photos show that
# Scale-free topology fit index as a function of the soft-thresholding power                                      
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power                                     
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


##Construct a co-expression matrix,according to best beta value
net = blockwiseModules(
  datExpr,
  power = sft$powerEstimate,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "AS-green-FPKM-TOM",
  verbose = 3
)
table(net$colors)
                                      

##Module visualization  
#Convert labels to colors for plotting                                      
mergedColors = labels2colors(net$colors)
table(mergedColors)
      
                                      
## Plot the dendrogram and the module colors underneath  
#assign all of the gene to their corresponding module 
# hclust for the genes.                                     
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
      
                                      
##相关性
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)                                     

                                      

#Clear sample size and number of genes                                      
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#First make a systematic clustering tree for the sample
datExpr_tree<-hclust(dist(datExpr), method = "average")
par(mar = c(0,5,2,0))
plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1, cex.main = 1,cex.lab=1)
#Add the corresponding color to the sample matrix of the previous construction
sample_colors <- numbers2colors(as.numeric(factor(datTraits$cellline)), 
                                colors = c("white","blue","red","green"),signed = FALSE)
#Constructing a system clustering tree and trait heat map of 10 samples
par(mar = c(1,4,3,1),cex=0.8)
plotDendroAndColors(datExpr_tree, sample_colors,
                    groupLabels = colnames(sample),
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")



##Relationship between modules and traits
design=model.matrix(~0+ datTraits$subtype)
colnames(design)=levels(datTraits$subtype)
moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0) #ME value matrix for modules of different colors (sample vs module)
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

                                      
##Specific genetic analysis of modules of interest traits                                    
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
## Calculate the Pearson correlation coefficient matrix for each module and gene
## MEs is the value of each module in each sample
## datExpr is the amount of expression of each gene in each sample
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

## Only continuous traits can only be calculated
## Here, the variable belonging to the Luminal phenotype is digitized with 0, 1
Luminal = as.data.frame(design[,3])
names(Luminal) = "Luminal"
geneTraitSignificance = as.data.frame(cor(datExpr, Luminal, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Luminal), sep="")
names(GSPvalue) = paste("p.GS.", names(Luminal), sep="")


module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Luminal",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

                                      
##Network visualization                                      
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
geneTree = net$dendrograms[[1]] 
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6)
plotTOM = dissTOM^7
diag(plotTOM) = NA
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")


# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
## Only continuous traits can only be calculated
## Here, the variable belonging to the Luminal phenotype is digitized with 0, 1
Luminal = as.data.frame(design[,3]);
names(Luminal) = "Luminal"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, Luminal))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
##Module clustering diagram
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
## Traits and module heat map
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


                                      
                                      
##Extract the gene name of the specified module                                      
# Select module
module = "brown" #Extract the corresponding color
# Select module probes
probes = colnames(datExpr) ## The probe in our example is the gene name.
inModule = (moduleColors==module)
modProbes = probes[inModule]


##Module export
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 6); 
# Select module
module = "brown";
# Select module probes
probes = colnames(datExpr) ## The probe in our example is the gene name.
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
## Also extract the gene name of the specified module
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
## Gene relationship matrix

#The first is to export to VisANT
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0)


#Then export to cytoscape
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes, 
  nodeAttr = moduleColors[inModule]
)

##If the module contains too many genes, the network is too complex and can be filtered.
nTop = 30
IMConn = softConnectivity(datExpr[, modProbes])
top = (rank(-IMConn) <= nTop)
filter <- modTOM[top, top]

                                      
##导出模块基因
TOM = TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate); 
# Select module
module = "turquoise"
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module)##导出基因名
inModule = datExpr1[,moduleColors==module] ##导出相应模块的基因表达矩阵                                     
modProbes = probes[inModule] 
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

##画出各组基因表达与特征值的趋势
par(mfrow=c(2,1), mar=c(0,4.1,4,2.05))
plotMat(t(scale(datExpr1[,moduleColors==module]) ),
        nrgcols=30,rlabels=F,rcols=module,
        main=module, cex.main=2)
par(mar=c(2,2.3,0.5,0.8))
barplot(as.matrix(t(MEs$MEbrown)), col=module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
                                      ####t(MEs$MEbrown)导出brown模块的各组特征值，做转置
