from subprocess import *
import sys
import os
import argparse
from os.path import basename

class Alignment():
    def __init__(self):
        self.samtools = "/home/zxchen/anaconda3/bin/samtools"
        self.hisat = "/home/zxchen/anaconda3/bin/hisat2"

    def run(self, cmd, wkdir=None):
        sys.stderr.write("Running %s ...\n" % cmd)
        p = Popen(cmd, shell=True, cwd=wkdir)
        p.wait()
        return p.returncode

    def fqTosam(self, txt, ref, fq_1,fq_2,log,sam):
        cmd = "%s --known-splicesite-infile %s -x %s --rna-strandness RF -1 %s -2 %s 2> %s -S %s" % (self.hisat,txt,ref,fq_1,fq_2,log,sam)
        return cmd

    def samTobam(self,sam,bam):
        cmd = "%s view -S %s -b > %s" % (self.samtools,sam,bam)
        return cmd

    def bamSort(self,bam,sortbam):
        cmd = "%s sort -n -o %s %s" % (self.samtools,sortbam,bam)
        return cmd

    def sortTormdup(self,sortbam, rmdup_bam):
        cmd = "%s rmdup %s %s" % (self.samtools,sortbam, rmdup_bam)
        return cmd


class Count():
    def __init__(self):
        self.HTseq = "/home/zxchen/anaconda3/bin/htseq-count"

    def htseq(self,rmdup_bam,gtf,htseq):
        cmd = "%s -f bam %s -r name -s reverse -t exon -i gene_id %s > %s" % (self.HTseq,rmdup_bam,gtf,htseq)
        return cmd

    def run(self, cmd, wkdir=None):
        sys.stderr.write("Running %s ...\n" % cmd)
        p = Popen(cmd, shell=True, cwd=wkdir)
        p.wait()
        return p.returncode


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="RNA")
    parser.add_argument("-ref", "--ref", required=True, type=str, help="the reference")
    parser.add_argument("-f1", "--f1", required=True, type=str, help="the fastq_1")
    parser.add_argument("-f2", "--f2", required=True, type=str, help="the fastq_2")
    parser.add_argument("-gtf", "--gtf", required=True, type=str, help="the gtf file")
    parser.add_argument("-txt", "--txt", required=True, type=str, help="the reference txt")
    parser.add_argument("-log", "--log", required=True, type=str, help="the log file")
    parser.add_argument("-ot", "--ot", required=True, type=str, help="the output dir")
    Args = parser.parse_args()

    step = 1
    al = Alignment()
    dir = Args.ot
    filename = basename(Args.f1)
    samfile = os.path.join(dir, 'sam/') + filename + '.sam'

    if step == 1:
        isExists = os.path.exists(os.path.join(dir,'sam'))
        if not isExists:
            os.mkdir(os.path.join(dir,'sam'))
            print('Create the dir!')
        else:
            print('Already exist!')

        cmd_1 = al.fqTosam(txt=Args.txt, ref=Args.ref, fq_1=Args.f1,fq_2=Args.f2,log=Args.log, sam=samfile)
        al.run(cmd=cmd_1)
        step = step + 1


    if step == 2:
        isExists = os.path.exists(os.path.join(dir,'bam'))
        if not isExists:
            os.mkdir(os.path.join(dir,'bam'))
            print('Create the dir!')
        else:
            print('Already exist!')

        bamfile = str(samfile).replace('sam', 'bam')
        cmd_2 = al.samTobam(sam=samfile,bam=bamfile)
        al.run(cmd=cmd_2)
        step = step + 1

    if step == 3:
        isExists = os.path.exists(os.path.join(dir,'sortbam'))
        if not isExists:
            os.mkdir(os.path.join(dir,'sortbam'))
            print('Create the dir!')
        else:
            print('Already exist!')
        bamfile = str(samfile).replace('sam', 'bam')
        sortbamfile = str(samfile).replace('sam', 'sortbam')
        cmd_3 = al.bamSort(sortbam=sortbamfile,bam=bamfile)
        al.run(cmd=cmd_3)
        step = step + 1

    if step == 4:
        isExists = os.path.exists(os.path.join(dir,'rmbam'))
        if not isExists:
            os.mkdir(os.path.join(dir,'rmbam'))
            print('Create the dir!')
        else:
            print('Already exist!')
        sortbamfile = str(samfile).replace('sam', 'sortbam')
        rmdupbamfile = str(samfile).replace('sam', 'rmbam')
        cmd_4 = al.sortTormdup(sortbam=sortbamfile, rmdup_bam=rmdupbamfile)
        al.run(cmd=cmd_4)
        step = step + 1

    cu = Count()

    if step == 5:
        isExists = os.path.exists(os.path.join(dir,'htseq'))
        if not isExists:
            os.mkdir(os.path.join(dir,'htseq'))
            print('Create the dir!')
        else:
            print('Already exist!')
        rmdupbamfile = str(samfile).replace('sam', 'rmbam')
        htseqfile = str(samfile).replace('sam', 'htseq')
        cmd_5 = cu.htseq(rmdup_bam=rmdupbamfile,gtf=Args.gtf,htseq=htseqfile)
        cu.run(cmd=cmd_5)











