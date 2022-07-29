BiocManager::install("maSigPro",version = "3.10",dependencies=TRUE, INSTALL_opts = c('--no-lock'))
install.packages("maSigPro",version = "3.10")
install.packages("venn")
install.packages("admisc")

install.packages("Cairo",dependencies=TRUE, INSTALL_opts = c('--no-lock'))



library(maSigPro)
data(data.abiotic)
data(edesign.abiotic)

design <- make.design.matrix(edesign.abiotic, degree = 2)

fit <- p.vector(data.abiotic, design, Q = 0.05, MT.adjust = "BH", min.obs = 20)
tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
sigs <- get.siggenes(tstep, rsq = 0.6, vars = "groups")

suma2Venn(sigs$summary[, c(2:4)])
suma2Venn(sigs$summary[, c(1:4)])

sigs$sig.genes$SaltvsControl$g
see.genes(sigs$sig.genes$ColdvsControl, 
          show.fit = T, dis =design$dis, 
          cluster.method="hclust" ,
          cluster.data = 1, k = 9)


STMDE66 <- data.abiotic[rownames(data.abiotic)=="STMDE66", ]
PlotGroups (STMDE66, edesign = edesign.abiotic)

PlotGroups (STMDE66, edesign = edesign.abiotic, 
            show.fit = T, dis = design$dis, 
            groups.vector = design$groups.vector)



cluster4 = data.frame(gene = c('ENSDARG00000010977',
                               'ENSDARG00000025671',
                               'ENSDARG00000033462',
                               'ENSDARG00000036344',
                               'ENSDARG00000037539',
                               'ENSDARG00000041399',
                               'ENSDARG00000069970',
                               'ENSDARG00000075780',
                               'ENSDARG00000075949',
                               'ENSDARG00000076103',
                               'ENSDARG00000077047',
                               'ENSDARG00000087191',
                               'ENSDARG00000093443',
                               'ENSDARG00000094559',
                               'ENSDARG00000095744',
                               'ENSDARG00000099470',
                               'ENSDARG00000102091',
                               'ENSDARG00000102364',
                               'ENSDARG00000110304',
                               'ENSDARG00000111040',
                               'ENSDARG00000112308',
                               'ENSDARG00000113526',
                               'ENSDARG00000113889',
                               'ENSDARG00000114226',
                               'ENSDARG00000116762',
                               'ENSDARG00000117485'
))

library(DOSE)
library(org.Dr.eg.db)
#library(topGO)
library(clusterProfiler)
#library(pathview)
library(AnnotationHub)

BP <- enrichGO(cluster4$gene,
                  "org.Dr.eg.db",
                  keyType="ENSEMBL",
                  ont="BP"
)

barplot(BP_up,title = "Biological process")

