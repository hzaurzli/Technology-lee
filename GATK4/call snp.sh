#bwa mapping
/home/yxtu/softwara/bwa/bwa-0.7.12/bwa mem -t 20 -M -R "@RG\tID:test\tSM:test\tLB:WGS\tPL:Illumina" 
"/home/yxtu/Ref_RNASeq/genome/zebrafish/NCBI/GCF_000002035.6_GRCz11_genomic.fna" 
"/home/yxtu/gRNA/Vazyme_240_S3_L001_R1_001.fastq" "/home/yxtu/gRNA/Vazyme_240_S3_L001_R2_001.fastq" > Vazyme_240_S3_L001.sam

#transform sam to sort.bam(sort by coordinate)
java -jar /home/yxtu/softwara/picard.jar SortSam I=/home/rzli/gRNA/Vazyme_240_S3_L001.sam 
o=/home/rzli/gRNA/Vazyme_240_S3_L001.sorted.bam SORT_ORDER=coordinate

#MarkDuplicates
nohup java -jar /home/yxtu/softwara/picard.jar MarkDuplicates  I=/home/rzli/gRNA/Vazyme_240_S3_L001.sorted.bam 
O=/home/rzli/gRNA/Vazyme_240_S3_L001.sorted_mark_dup.bam M=/home/rzli/gRNA/Vazyme_240_S3_L001.sorted_mark_dup.txt &

##three files in same dirction
#creat dic
nohup java -jar /home/yxtu/softwara/picard.jar CreateSequenceDictionary REFERENCE=/home/rzli/gRNA/gatk/GCF_000002035.6_GRCz11_genomic.fna 
OUTPUT=/home/rzli/gRNA/GCF_000002035.6_GRCz11_genomic.dict &

#index
samtools faidx GCF_000002035.6_GRCz11_genomic.fna

##gatk check var
#index dic fa in same dirction
nohup java -jar /home/dywang/software/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar HaplotypeCaller 
-R "/home/rzli/gRNA/gatk/GCF_000002035.6_GRCz11_genomic.fna" -I /home/rzli/gRNA/Vazyme_240_S3_L001.sorted_mark_dup_add.bam 
-O /home/rzli/gRNA/Vazyme_240_S3_L001.sorted_mark_dup_add.g.vcf.gz -ERC GVCF &
