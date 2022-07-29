#cuffquant
~/biosoft/cufflinks-2.2.1.Linux_x86_64/cuffquant -o ~/swxxxjz/mapping/tot/sex.gtf/1.0_sex/1.0_sex_quant -p 8 -u ~/swxxxjz/Danio_gff3 ~/swxxxjz/mapping/tot/sex.gtf/1.0_sex/1.0_sex.bam

#cuffnorm
~/biosoft/cufflinks-2.2.1.Linux_x86_64/cuffnorm -o ~/swxxxjz/mapping/tot/sex.gtf/0.5_sex/cuffnorm_out -p 8 -L 0.5_sex,1.0_sex ~/swxxxjz/Danio_gff3 ~/swxxxjz/mapping/tot/sex.gtf/0.5_sex/0.5_sex.quant/abundances.cxb ~/swxxxjz/mapping/tot/sex.gtf/1.0_sex/1.0_sex_quant/abundances.cxb


#####循环
 for i in $(ls *bam ~/swxxxjz/mapping/tot/sex.gtf/*/*bam);do echo ~/biosoft/cufflinks-2.2.1.Linux_x86_64/cuffquant -o ${i%.bam}.quant -p 8 -u ~/swxxxjz/Danio_gff3 $i;done


######sex
~/biosoft/cufflinks-2.2.1.Linux_x86_64/cuffnorm -o ~/swxxxjz/bam/tot/sex.gtf/cuffnorm_out -p 8 -L 0.1_sex,0.2_sex,0.3_sex,0.4_sex,0.5_sex,0.6_sex,0.7_sex,0.8_sex,0.9_sex,1.0_sex ~/swxxxjz/Danio_gff3 ~/swxxxjz/bam/tot/sex.gtf/0.1_sex/0.1_sex.quant/abundances.cxb ~/swxxxjz/bam/tot/sex.gtf/0.2_sex/0.2_sex.quant/abundances.cxb ~/swxxxjz/bam/tot/sex.gtf/0.3_sex/0.3_sex.quant/abundances.cxb ~/swxxxjz/bam/tot/sex.gtf/0.4_sex/0.4_sex.quant/abundances.cxb ~/swxxxjz/bam/tot/sex.gtf/0.5_sex/0.5_sex.quant/abundances.cxb ~/swxxxjz/bam/tot/sex.gtf/0.6_sex/0.6_sex.quant/abundances.cxb ~/swxxxjz/bam/tot/sex.gtf/0.7_sex/0.7_sex.quant/abundances.cxb ~/swxxxjz/bam/tot/sex.gtf/0.8_sex/0.8_sex.quant/abundances.cxb ~/swxxxjz/bam/tot/sex.gtf/0.9_sex/0.9_sex.quant/abundances.cxb ~/swxxxjz/bam/tot/sex.gtf/1.0_sex/1.0_sex.quant/abundances.cxb

######m
~/biosoft/cufflinks-2.2.1.Linux_x86_64/cuffnorm -o ~/swxxxjz/bam/tot/m.gtf/cuffnorm_out -p 8 -L 0.1_m,0.2_m,0.3_m,0.4_m,0.5_m,0.6_m,0.7_m,0.8_m,0.9_m,1.0_m ~/swxxxjz/Danio_gff3 ~/swxxxjz/bam/tot/m.gtf/0.1_m/0.1_m.quant/abundances.cxb ~/swxxxjz/bam/tot/m.gtf/0.2_m/0.2_m.quant/abundances.cxb ~/swxxxjz/bam/tot/m.gtf/0.3_m/0.3_m.quant/abundances.cxb ~/swxxxjz/bam/tot/m.gtf/0.4_m/0.4_m.quant/abundances.cxb ~/swxxxjz/bam/tot/m.gtf/0.5_m/0.5_m.quant/abundances.cxb ~/swxxxjz/bam/tot/m.gtf/0.6_m/0.6_m.quant/abundances.cxb ~/swxxxjz/bam/tot/m.gtf/0.7_m/0.7_m.quant/abundances.cxb ~/swxxxjz/bam/tot/m.gtf/0.8_m/0.8_m.quant/abundances.cxb ~/swxxxjz/bam/tot/m.gtf/0.9_m/0.9_m.quant/abundances.cxb ~/swxxxjz/bam/tot/m.gtf/1.0_m/1.0_m.quant/abundances.cxb


######6
~/biosoft/cufflinks-2.2.1.Linux_x86_64/cuffnorm -o ~/swxxxjz/bam/tot/6.gtf/cuffnorm_out -p 8 -L 0.1_6,0.2_6,0.3_6,0.4_6,0.5_6,0.6_6,0.7_6,0.8_6,0.9_6,1.0_6 ~/swxxxjz/Danio_gff3 ~/swxxxjz/bam/tot/6.gtf/0.1_6/0.1_6.quant/abundances.cxb ~/swxxxjz/bam/tot/6.gtf/0.2_6/0.2_6.quant/abundances.cxb ~/swxxxjz/bam/tot/6.gtf/0.3_6/0.3_6.quant/abundances.cxb ~/swxxxjz/bam/tot/6.gtf/0.4_6/0.4_6.quant/abundances.cxb ~/swxxxjz/bam/tot/6.gtf/0.5_6/0.5_6.quant/abundances.cxb ~/swxxxjz/bam/tot/6.gtf/0.6_6/0.6_6.quant/abundances.cxb ~/swxxxjz/bam/tot/6.gtf/0.7_6/0.7_6.quant/abundances.cxb ~/swxxxjz/bam/tot/6.gtf/0.8_6/0.8_6.quant/abundances.cxb ~/swxxxjz/bam/tot/6.gtf/0.9_6/0.9_6.quant/abundances.cxb ~/swxxxjz/bam/tot/6.gtf/1.0_6/1.0_6.quant/abundances.cxb


######51
~/biosoft/cufflinks-2.2.1.Linux_x86_64/cuffnorm -o ~/swxxxjz/bam/tot/51.gtf/cuffnorm_out -p 8 -L 0.1_51,0.2_51,0.3_51,0.4_51,0.5_51,0.6_51,0.7_51,0.8_51,0.9_51,1.0_51 ~/swxxxjz/Danio_gff3 ~/swxxxjz/bam/tot/51.gtf/0.1_51/0.1_51.quant/abundances.cxb ~/swxxxjz/bam/tot/51.gtf/0.2_51/0.2_51.quant/abundances.cxb ~/swxxxjz/bam/tot/51.gtf/0.3_51/0.3_51.quant/abundances.cxb ~/swxxxjz/bam/tot/51.gtf/0.4_51/0.4_51.quant/abundances.cxb ~/swxxxjz/bam/tot/51.gtf/0.5_51/0.5_51.quant/abundances.cxb ~/swxxxjz/bam/tot/51.gtf/0.6_51/0.6_51.quant/abundances.cxb ~/swxxxjz/bam/tot/51.gtf/0.7_51/0.7_51.quant/abundances.cxb ~/swxxxjz/bam/tot/51.gtf/0.8_51/0.8_51.quant/abundances.cxb ~/swxxxjz/bam/tot/51.gtf/0.9_51/0.9_51.quant/abundances.cxb ~/swxxxjz/bam/tot/51.gtf/1.0_51/1.0_51.quant/abundances.cxb

#############################
#############################
#seqtk
~/biosoft/seqtk-1.2/seqtk sample -s100 51_clean_r2.fq 28882130 > ~/swxxxjz/80/51.tes_r2.fq

#bam随机提取50%bam
samtools view -s 0.5 -bo out.bam in.bam

#stringtie转录本定量
stringtie out.bam -o test.gtf -p 20 -G /home/yxtu/Ref_RNASeq/annotation/zebrafish/Ensembl_96/Danio_rerio.GRCz11.96.gtf

#merge转录本
stringtie --merge -p 8 -G /home/yxtu/Ref_RNASeq/annotation/zebrafish/Ensembl_96/Danio_rerio.GRCz11.96.gtf -o merge_m.gtf ~/swxxxjz/mapping/tot/m.gtf/mergelist.txt

#转录本矫正
gffcompare -r /home/yxtu/Ref_RNASeq/annotation/zebrafish/Ensembl_96/Danio_rerio.GRCz11.96.gtf -G -o merged merge_m.gtf

#重新定量
stringtie -e -B -p 8 -G merge_m.gtf -o 0.1_m.test_a.gtf 0.1_m.bam

stringtie -e -B -p 8 -G test.gtf -o test_output.gtf input.bam –A gene_out

for i in *m;do echo stringtie -e -B -p 8 -G ~/swxxxjz/mapping/tot/m.gtf/merge_m.gtf -o ~/swxxxjz/mapping/tot/m.gtf/$i/${i}.test_a.gtf ~/swxxxjz/mapping/tot/m.gtf/*/${i}.bam ;done

#循环stringtie转录本定量
for i in $(ls *bam ~/swxxxjz/bam/tot/m.gtf/*/*bam);do echo stringtie $i -o ${i%.bam}.gtf -p 8 -G /home/yxtu/Ref_RNASeq/annotation/zebrafish/Ensembl_96/Danio_rerio.GRCz11.96.gtf;done

stringtie --merge -p 8 -G /home/yxtu/Ref_RNASeq/annotation/zebrafish/Ensembl_96/Danio_rerio.GRCz11.96.gtf -o merge_m.gtf ~/swxxxjz/bam/tot/sex.gtf/stringmerge/mergelist.txt



for i in $(ls *m.bam ~/swxxxjz/bam/tot/m.gtf/*/*m.gtf);do mv $i ~/swxxxjz/bam/tot/m.gtf/stringmerge/ ;done
