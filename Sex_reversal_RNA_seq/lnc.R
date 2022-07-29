gtf1 <- rtracklayer::import('/home/ug0007/rzli/zebrafish.gtf')
gtf_df <- as.data.frame(gtf1)
test <- gtf_df[1:20,]
View(test)


gtf_se = gtf_df[,c(7,10,12,14)]
gtf_gene = subset(gtf_se,gtf_se[,1] == "gene")


table(gtf_gene$gene_biotype)

gtf_lnc = subset(  
  gtf_gene,gtf_gene[,4] == "lincRNA" | gtf_gene[,4] == "pseudogene" | 
    gtf_gene[,4] == "IG_pseudogene" | gtf_gene[,4] == "IG_J_pseudogene" |
    gtf_gene[,4] == "IG_V_pseudogene" | gtf_gene[,4] == "processed_pseudogene" |
    gtf_gene[,4] == "TR_V_pseudogene" |  gtf_gene[,4] == "polymorphic_pseudogene" | 
    gtf_gene[,4] == "processed_transcript" | gtf_gene[,4] == "transcribed_unprocessed_pseudogene")

gtf_coding = subset(gtf_gene,gtf_gene[,4] == "protein_coding" )
