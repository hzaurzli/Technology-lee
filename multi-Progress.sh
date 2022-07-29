#!/bin/bash
cat /home/rzli/swxxxjz/bam/tot/13.gtf/select | while read line
do
 {
  $line
 }&
done
wait


##The file select
samtools view -s 0.02 -bo /home/rzli/swxxxjz/bam/tot/13.gtf/0.02_13/0.02_13.bam /home/rzli/swxxxjz/bam/tot/13.gtf/1.0_13.bam
samtools view -s 0.04 -bo /home/rzli/swxxxjz/bam/tot/13.gtf/0.04_13/0.04_13.bam /home/rzli/swxxxjz/bam/tot/13.gtf/1.0_13.bam
samtools view -s 0.06 -bo /home/rzli/swxxxjz/bam/tot/13.gtf/0.06_13/0.06_13.bam /home/rzli/swxxxjz/bam/tot/13.gtf/1.0_13.bam
samtools view -s 0.08 -bo /home/rzli/swxxxjz/bam/tot/13.gtf/0.08_13/0.08_13.bam /home/rzli/swxxxjz/bam/tot/13.gtf/1.0_13.bam
samtools view -s 0.1 -bo /home/rzli/swxxxjz/bam/tot/13.gtf/0.1_13/0.1_13.bam /home/rzli/swxxxjz/bam/tot/13.gtf/1.0_13.bam
samtools view -s 0.2 -bo /home/rzli/swxxxjz/bam/tot/13.gtf/0.2_13/0.2_13.bam /home/rzli/swxxxjz/bam/tot/13.gtf/1.0_13.bam
samtools view -s 0.3 -bo /home/rzli/swxxxjz/bam/tot/13.gtf/0.3_13/0.3_13.bam /home/rzli/swxxxjz/bam/tot/13.gtf/1.0_13.bam
samtools view -s 0.4 -bo /home/rzli/swxxxjz/bam/tot/13.gtf/0.4_13/0.4_13.bam /home/rzli/swxxxjz/bam/tot/13.gtf/1.0_13.bam
samtools view -s 0.5 -bo /home/rzli/swxxxjz/bam/tot/13.gtf/0.5_13/0.5_13.bam /home/rzli/swxxxjz/bam/tot/13.gtf/1.0_13.bam
samtools view -s 0.6 -bo /home/rzli/swxxxjz/bam/tot/13.gtf/0.6_13/0.6_13.bam /home/rzli/swxxxjz/bam/tot/13.gtf/1.0_13.bam
samtools view -s 0.7 -bo /home/rzli/swxxxjz/bam/tot/13.gtf/0.7_13/0.7_13.bam /home/rzli/swxxxjz/bam/tot/13.gtf/1.0_13.bam
samtools view -s 0.8 -bo /home/rzli/swxxxjz/bam/tot/13.gtf/0.8_13/0.8_13.bam /home/rzli/swxxxjz/bam/tot/13.gtf/1.0_13.bam
samtools view -s 0.9 -bo /home/rzli/swxxxjz/bam/tot/13.gtf/0.9_13/0.9_13.bam /home/rzli/swxxxjz/bam/tot/13.gtf/1.0_13.bam
samtools view -s 1.0 -bo /home/rzli/swxxxjz/bam/tot/13.gtf/1.0_13/1.0_13.bam /home/rzli/swxxxjz/bam/tot/13.gtf/1.0_13.bam


