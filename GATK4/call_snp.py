from subprocess import *
import sys
import os
import argparse


class callVars():
    def __init__(self):
        self.samtools = "/home/zxchen/anaconda3/bin/samtools"
        self.gatk = "/home/rzli/biosoft/gatk-4.1.5.0/gatk"
        self.bwa = "/home/yxtu/software/bwa/bwa-0.7.12/bwa"

    def run(self, cmd, wkdir=None):
        sys.stderr.write("Running %s ...\n" % cmd)
        p = Popen(cmd, shell=True, cwd=wkdir)
        p.wait()
        return p.returncode

    def bwa_index(self, ref=None):
        cmd = "%s index %s" % (self.bwa, ref)
        return cmd

    def samtools_index(self, ref=None):
        cmd = "%s faidx %s" % (self.samtools, ref)
        return cmd

    def gatk_CreateSequenceDictionary(self, ref=None, ref_dic=None):
        cmd = "%s CreateSequenceDictionary -R %s -O %s" % (self.gatk, ref, ref_dic)
        return cmd

    def bwa_alignment(self, ref=None, test_1_fq=None, test_2_fq=None, sam=None):
        cmd = "%s mem -t 4 -R %s %s %s >%s" % (self.bwa, ref, test_1_fq, test_2_fq,sam)
        return cmd

    def samtobam(self, sam=None, bam=None):
        cmd = "%s view -b -S %s > %s" % (self.samtools, sam, bam)
        return cmd

    def samtools_sort(self, sort_bam=None, bam=None):
        cmd = "%s sort -@ 3 -o %s %s" % (self.samtools, sort_bam, bam)
        return cmd

    def gatk_AddOrReplaceReadGroups(self, sort_bam=None, add_bam=None):
        cmd = "%s AddOrReplaceReadGroups -I %s -O %s -LB library1 -PL illumina -PU pl1 -SM name" % (self.gatk, sort_bam, add_bam)
        return cmd

    def gatk_MarkDuplicates(self, add_bam=None, remark_bam=None ,M=None):
        cmd = "%s MarkDuplicates -I %s -O %s -M %s" % (self.gatk, add_bam, remark_bam, M)
        return cmd

    def gatk_call_snp(self, remark_bam=None, vcf=None, ref=None):
        cmd = "%s --java-options -Xmx4G HaplotypeCaller -I %s -O %s -R %s" % (self.gatk, remark_bam, vcf, ref)
        return cmd

    def gatk_SelectVariants(self, vcf=None, indel_vcf=None, Type=None):
        cmd = "%s SelectVariants -V %s -O %s --select-type-to-include %s" % (self.gatk,vcf ,indel_vcf , Type)
        return cmd

class annotation():
    def __init__(self):
        self.gtfToGenePred = "/home/rzli/biosoft/gtfToGenePred"
        self.retrieve_seq = "/home/rzli/biosoft/annovar/retrieve_seq_from_fasta.pl"
        self.convert2annovar = "/home/rzli/biosoft/annovar/convert2annovar.pl"
        self.annotate_variation = "/home/rzli/biosoft/annovar/annotate_variation.pl"

    def run(self, cmd, wkdir=None):
        sys.stderr.write("Running %s ...\n" % cmd)
        p = Popen(cmd, shell=True, cwd=wkdir)
        p.wait()
        return p.returncode

    def gtfToGene(self, gtf=None, ref_txt=None):
        cmd = "%s -genePredExt %s %s" % (self.gtfToGenePred ,gtf, ref_txt)
        return cmd

    def retrieve(self,fa=None ,ref_txt=None ,rna_fa=None):
        cmd = "%s --format refGene --seqfile %s %s --out %s" % (self.retrieve_seq, fa, ref_txt, rna_fa)
        return cmd

    def convert(self, vcf=None, O=None):
        cmd = "%s -format vcf4old %s >%s" % (self.convert2annovar, vcf, O)
        return cmd

    def annotate_variation(self,type=None, O=None, input=None):
        cmd = "%s -geneanno --neargene 2000 -buildver %s -dbtype refGene -outfile %s -exonsort %s" % (self.annotate_variation, type, O, input)
        return cmd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="SNP")
    parser.add_argument("-ref", "--ref", required=True, type=str, help="the reference")
    parser.add_argument("-f1", "--f1", required=True, type=str, help="the fastq1")
    parser.add_argument("-f2", "--f2", required=True, type=str, help="the fastq2")
    parser.add_argument("-gtf", "--gtf", required=True, type=str, help="the gtf file")

    Args = parser.parse_args()
    cv = callVars()
    #os.mkdir("/home/rzli/vars_new")

    # step 1
    # os.chdir("/home/rzli/vars_new/")
    # os.mkdir("/home/rzli/vars_new/ref")
    cmd_1 = cv.bwa_index(ref=Args.ref)
    cv.run(cmd=cmd_1)
    #step 2
    os.chdir("/home/rzli/vars_new/")
    os.mkdir("/home/rzli/vars_new/bam")
    cmd_2 = cv.samtools_index(ref=Args.ref)
    cv.run(cmd=cmd_2)
    #step 3
    os.chdir("/home/rzli/vars_new/")
    filepath = "/home/rzli/vars_new/ref.dict"
    cmd_3 = cv.gatk_CreateSequenceDictionary(ref=Args.ref ,ref_dic=filepath)
    cv.run(cmd=cmd_3)
    #step 4
    os.chdir("/home/rzli/vars_new/bam/")
    bampath = "/home/rzli/vars_new/bam/test.bam"
    cmd_4 = cv.bwa_alignment(ref=Args.ref, test_1_fq=Args.f1, test_2_fq=Args.f2, sam=sampath)
    cv.run(cmd=cmd_4)
    #step 5
    os.chdir("/home/rzli/vars_new/bam/")
    sortpath = "/home/rzli/vars_new/bam/test_sort.bam"
    cmd_5 = cv.samtools_sort(sort_bam=sortpath, bam=bampath)
    cv.run(cmd=cmd_5)
    #step 6
    os.chdir("/home/rzli/vars_new/bam/")
    add_bam = "/home/rzli/vars_new/bam/add.bam"
    cmd_6 = cv.gatk_AddOrReplaceReadGroups(sort_bam=sortpath, add_bam=add_bam)
    cv.run(cmd=cmd_6)
    #step 7
    os.chdir("/home/rzli/vars_new/bam/")
    remark_bam = "/home/rzli/vars_new/bam/remark.bam"
    Metrics = "/home/rzli/vars_new/bam/markdup_metrics.txt"
    cmd_7 = cv.gatk_MarkDuplicates(add_bam=add_bam, remark_bam=remark_bam ,M=Metrics)
    cv.run(cmd=cmd_7)
    #step 8
    os.chdir("/home/rzli/vars_new/")
    os.mkdir("/home/rzli/vars_new/vcf")
    vcf = "/home/rzli/vars_new/vcf/test.vcf"
    cmd_8 = cv.gatk_call_snp(remark_bam=remark_bam, vcf=vcf, ref=Args.r)
    cv.run(cmd=cmd_8)
    #step 9
    os.chdir("/home/rzli/vars_new/vcf")
    indel_vcf = "/home/rzli/vars_new/vcf/test_indel.vcf"
    Type = "INDEL"
    cmd_9 = cv.gatk_SelectVariants(vcf=vcf, indel_vcf=indel_vcf, Type=Type)
    cv.run(cmd=cmd_9)

    #step 10
    ann = annotation()
    os.chdir("/home/rzli/vars_new/")
    os.mkdir("/home/rzli/vars_new/annotation")
    gtf = Args.gtf
    ref_txt = "/home/rzli/vars_new/annotation/ref_txt"
    cmd_10 = ann.gtfToGene(gtf=gtf, ref_txt=ref_txt)
    ann.run(cmd=cmd_10)

    # step 11
    os.chdir("/home/rzli/vars_new/annotation/")
    fa = Args.ref
    rna_fa = "/home/rzli/vars_new/annotation/genome_refGeneMrna.fa"
    cmd_11 = ann.retrieve(fa=fa ,ref_txt=ref_txt ,rna_fa=rna_fa)
    ann.run(cmd=cmd_11)

    #step 12
    os.chdir("/home/rzli/vars_new/annotation/")
    out = "/home/rzli/vars_new/annotation/annotation_input"
    cmd_12 = ann.convert(vcf=indel_vcf, O=out)
    ann.run(cmd=cmd_12)





