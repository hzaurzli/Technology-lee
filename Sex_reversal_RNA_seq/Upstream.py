import os
import subprocess as sub

def qsub(pbsfile):
    cmd = "qsub " + pbsfile
    sub.call(cmd, shell=True)

samples = []
for i in os.listdir("/home/rzli/fish/1"):
    i = i.split("_")
    sample = i[0] + "_" + i[1]
    samples.append(sample)
    samples = list(set(samples))


for i in samples:
    with open("/home/rzli/fish/zreb/codes/" + i + ".sh", "w") as w:
        line = "#!/bin/bash" + "\n"
        line = line + "hisat2 --known-splicesite-infile /home/yxtu/Ref_RNASeq/annotation/zebrafish/Ensembl_96/zebrafish_96.txt -x /home/yxtu/Ref_RNASeq/genome/zebrafish/genome_96/genome --rna-strandness RF "
        line = line + "-1 /home/rzli/fish/1/" + i + "_1.clean.fq.gz "
        line = line + "-2 /home/rzli/fish/2/" + i + "_2.clean.fq.gz "
        line = line + "2>/home/rzli/fish/zreb/log/" + i + ".log" + " | samtools view - -Sb | samtools sort - -T "
        line = line + "/home/rzli/fish/zreb/mapping/" + i + " -o " + "/home/rzli/fish/zreb/mapping/" + i + ".bam"
        w.write(line)

for i in samples:
    qsub("/home/rzli/fish/zreb/codes/" + i + ".sh")
    
    
samples = []
for i in os.listdir("/home/rzli/fish/zreb/mapping"):
    i = i.split(".")
    sample = i[0]
    samples.append(sample)
    samples = list(set(samples))

for i in samples:
    with open("/home/rzli/fish/zreb/htcodes/" + i + ".sh", "w") as w:
        line = "#!/bin/bash" + "\n"
        line = line + "htseq-count -f bam /home/rzli/fish/zreb/mapping/" + i + ".bam " + "-r name -s reverse -t exon -i gene_id /home/yxtu/Ref_RNASeq/annotation/zebrafish/Ensembl_96/Danio_rerio.GRCz11.96.gtf > /home/rzli/fish/zreb/htseq/" + i +".htseq"
        w.write(line)

for i in samples:
    qsub("/home/rzli/fish/zreb/htcodes/" + i + ".sh")
