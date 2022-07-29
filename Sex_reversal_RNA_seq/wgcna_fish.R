####################

expro = cbind(data_64_my,data_32,data_lav)
expro = apply(expro,2,function(x) {x/sum(x)*1000000})

expro = cbind(res_32,res_64,res_lav)

expro = expro[which(row.names(expro) %in% exp$select_genes),]


samples = data.frame(
  DMSO_32_F = c(rep(1,2),rep(0,27)),
  DMSO_32_M = c(rep(0,2),rep(1,2),rep(0,25)),
  EM_32_F = c(rep(0,4),rep(1,2),rep(0,23)),
  EM_32_M = c(rep(0,6),rep(1,2),rep(0,21)),
  DMSO_64_F = c(rep(0,8),rep(1,4),rep(0,17)),
  DMSO_64_M = c(rep(0,12),rep(1,2),rep(0,15)),
  EM_64_F = c(rep(0,14),rep(1,3),rep(0,12)),
  EM_64_M = c(rep(0,17),rep(1,2),rep(0,10)),
  DMSO_32_lav = c(rep(0,19),rep(1,3),rep(0,7)),
  EM_32_lav = c(rep(0,22),rep(1,3),rep(0,4)),
  DMSO_64_lav = c(rep(0,25),rep(1,2),rep(0,2)),
  EM_64_lav = c(rep(0,27),rep(1,2))
)


expro = expro[,-c(9:19)]
samples = data.frame(
  DMSO_32_F = c(rep(1,2),rep(0,16)),
  DMSO_32_M = c(rep(0,2),rep(1,2),rep(0,14)),
  EM_32_F = c(rep(0,4),rep(1,2),rep(0,12)),
  EM_32_M = c(rep(0,6),rep(1,2),rep(0,10)),
  DMSO_32_lav = c(rep(0,8),rep(1,3),rep(0,7)),
  EM_32_lav = c(rep(0,11),rep(1,3),rep(0,4)),
  DMSO_64_lav = c(rep(0,14),rep(1,2),rep(0,2)),
  EM_64_lav = c(rep(0,16),rep(1,2))
)


row.names(samples) = colnames(expro)

#################################
datExpr=as.data.frame(t(expro));
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

powers = c(c(1:10), seq(from = 8, to=30, by=2))
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
net = blockwiseModules(datExpr, power = 28, deepSplit = 4,maxBlockSize =7455,
                       TOMType = "unsigned",minModuleSize = 100,
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


##表型与模块相关性##
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

