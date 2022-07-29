##导入数据##
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

samples=read.csv('/home/ug0007/rzli/bio/Sam_info.txt',sep = '\t',row.names = 1)
expro=read.csv('/home/ug0007/rzli/bio/ExpData.txt',sep = '\t',row.names = 1)
dim(expro)

##筛选方差前25%的基因##
m.vars=apply(expro,1,var)
expro.upper=expro[which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[4]),]
dim(expro.upper)
datExpr=as.data.frame(t(expro.upper));
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

##样本聚类检查离群值##
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK
sampleTree = hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers"
     , sub="", xlab="")
save(datExpr, file = "FPKM-01-dataInput.RData")

##没有离群值，不作处理##

##软阈值筛选##
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##一步法网络构建：One-step network construction and module detection##
net = blockwiseModules(datExpr, power = 14, maxBlockSize = 6000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "AS-green-FPKM-TOM",
                       verbose = 3)
table(net$colors)

##绘画结果展示##
# open a graphics window
#sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

##结果保存
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "AS-green-FPKM-02-networkConstruction-auto.RData")


##表型与模块相关性##
trainDt= allTraits[Rat, -1]
rownames(trainDt) = samples[traitRows, 1]
table(rownames(trainDt) == rownames(datExpr))

moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01, 
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))