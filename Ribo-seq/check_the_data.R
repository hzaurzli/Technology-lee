if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19.masked")


source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")


library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


install.packages("Biostrings")
install.packages("BiocGenerics")
install.packages("parallel")
install.packages("S4Vectors")
install.packages("stats4")
install.packages("IRanges")
install.packages("XVector")

install.packages("GenomicFeatures")
install.packages("GenomeInfoDb")
install.packages("GenomicRanges")
install.packages("AnnotationDbi")
install.packages("Biobase")

install.packages("rtracklayer")


data(ctrlGAlignments)
aln <- ctrlGAlignments
matchLenDistr <- histMatchLength(aln, 0)
matchLenDistr[[2]]

alnGRanges <- readsToStartOrEnd(aln, what="start")
#txdb object with annotations
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
oneBinRanges <- aroundPromoter(txdb, alnGRanges, percBestExpressed=0.001)
#the coverage in the TSS flanking region for the reads with match sizes 29:31
listPromoterCov <-
  readStartCov(
    alnGRanges,
    oneBinRanges,
    matchSize=c(29:31),
    fixedInterval=c(-20, 20),
    renameChr="aroundTSS",
    charPerc="perc"
  )
plotSummarizedCov(listPromoterCov)


alnGRanges <- alnGRanges[which(!is.na(match(alnGRanges$score,30:33)))]
#get all CDSs by transcript
cds <- cdsBy(txdb, by="tx", use.names=TRUE)
#get all exons by transcript
exonGRanges <- exonsBy(txdb, by="tx", use.names=TRUE)
#get the per transcript relative position of start and end codons
cdsPosTransc <- orfRelativePos(cds, exonGRanges)
#compute the counts on the different features
#after applying the specified shift value on the read start along the transcript
countsDataCtrl1 <-
  countShiftReads(
    exonGRanges=exonGRanges[names(cdsPosTransc)],
    cdsPosTransc=cdsPosTransc,
    alnGRanges=alnGRanges,
    shiftValue=-14
  )
head(countsDataCtrl1[[1]])
listCountsPlots <- countsPlot(
  list(countsDataCtrl1[[1]]),
  grep("_counts$", colnames(countsDataCtrl1[[1]])),
  1
)
listCountsPlots


data(codonIndexCovCtrl)
listReadsCodon <- countsDataCtrl1[[2]]
#get the names of the expressed ORFs grouped by transcript
cds <- cdsBy(txdb, use.names=TRUE)
orfCoord <- cds[names(cds) %in% names(listReadsCodon)]
#chromosome names should correspond between the BAM,
#the annotations, and the genome
genomeSeq <- BSgenome.Hsapiens.UCSC.hg19
codonData <- codonInfo(listReadsCodon, genomeSeq, orfCoord)

codonUsage <- codonData[[1]]
codonCovMatrix <- codonData[[2]]
#keep only genes with a minimum number of reads
nbrReadsGene <- apply(codonCovMatrix, 1, sum)
ixExpGenes <- which(nbrReadsGene >= 50)
codonCovMatrix <- codonCovMatrix[ixExpGenes, ]
#get the PCA on the codon coverage
codonCovMatrixTransp <- t(codonCovMatrix)
rownames(codonCovMatrixTransp) <- colnames(codonCovMatrix)
colnames(codonCovMatrixTransp) <- rownames(codonCovMatrix)
listPCACodonCoverage <- codonPCA(codonCovMatrixTransp, "codonCoverage")
listPCACodonCoverage[[2]]

