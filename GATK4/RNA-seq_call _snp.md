```
使用BWA软件(也是此文作者使用)对基因组进行建立索引。
bwa index coconut_dwarf.fasta

 
 使用bwa进行比对
bwa mem ref.fa sample_1.fq sample_2.fq -R '@RG\tID:sample\tLB:sample\tSM:sample\tPL:ILLUMINA' \
2>sample_map.log | samtools sort -@ 20 -O bam -o sample.sorted.bam 1>sample_sort.log 2>&1

#call snp
java -jar ~/packages/picard-tools-1.124/picard.jar MarkDuplicates INPUT=D1801357-GQ_S1_L008.sorted.bam OUTPUT=D1801357-GQ_S1_L008_sorted__dupl2.bam METRICS_FILE=D1801357-GQ_S1_L008_sorted_metric_2.txt#标记重复

~/packages/samtools-1.9/samtools mpileup -ug -f ~/public/data/Genome/tall/coconut_tall.fasta dwarf.sorted.bam | ~/packages/bcftools/bcftools call -mv --output-type v > dwarf.raw.vcf




```

