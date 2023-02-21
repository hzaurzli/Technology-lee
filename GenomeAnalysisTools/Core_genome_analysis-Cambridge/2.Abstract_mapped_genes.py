# This is the script for abstract mapped genes from blast files.(the extension is .sc19, because it was used for anaylsis the Streptococcus suis, and the SC19 isoaltes was chosen as the reference genome. You can change it in the script...do not forget the blast output formate is 6).
# Just put the script with all your blast files as well as all your genome files (the extension is .fna) under a same folder, and run the script, you will got the mapped-genes-abstracted file for each genome (the output file are named with adding "_sc19_blast.fasta", you can change it). The cutoff values for determining the mapped genes are 80% alignment length coverage and 80% identity. (Notic that, in order to calculat the alignment length coverage, your reference-genome-ffn-file blast database should be modified that adding the gene length to the ID of each genes behind a colon ":", if you do not know how to do it, there is a script named Add_seqlen_to_title.py will help you.)
# This script will also output a file named "data_results.out" that list the mapped genes for each genome. Maybe you need to handle the output file by yourself to creat the core genome genes list and the isolates list for further analysis, I do not crate a script for this step, it won't be hard...hopefully.
# This is just a home used script, not be polished, very personal style and no debugging, maybe not suitable for public using. Writtern by Geng ZOU,on 2022-11-02, if you have any questions, contact me. zougeng19900918@126.com


#-*-coding:utf-8-*-
import re
import sys, getopt
import operator
import os
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def blast_filter():
    ARG_location_dict = {}
    rootdir1 = "./"
    for parent,dirnames,filenames in os.walk(rootdir1):
        for filename in filenames:
            cc = filename
            dd = cc.split(".")[0]
            ff = cc.split(".")[-1]
            if ff == "sc19" : # the extension of the blast files
                blast_info = open(filename, "r")
                for line in blast_info:
                    line_info = line.strip().split("\t")
                    Contig_ID_info = line_info[0]
                    ARG_ID = line_info[1]
                    ARG_len = line_info[1].split(":")[-1]
                    identical_percent = line_info[2]
                    align_length = int(line_info[3])
                    contig_start = int(line_info[6])
                    contig_end = int(line_info[7])
                    ARG_start = int(line_info[8])
                    ARG_end = int(line_info[9])
                    Score = float(line_info[11])
                    F_R = "F"
                    if ARG_end < ARG_start:
                        F_R = "R"
                    align_percent = '%.2f' %(float(align_length)/float(ARG_len)*100)
                    key_use = (dd, Contig_ID_info)
                    if float(align_percent) > 50.00:
                        ARG_location_dict.setdefault(key_use,[]).append((ARG_ID, contig_start, contig_end, Score, identical_percent, align_percent, Contig_ID_info, F_R))
                blast_info.close()
    return ARG_location_dict

def ARG_filter(ARG_location_dict):
    ARG_location_filter_dict = {}
    ARG_list_get = []
    for item in ARG_location_dict.items():
        key_data = item[0]
        ARG_list = item[1]
        ARG_list.sort(key = operator.itemgetter(1))
        start_initial = 0
        ii_keep = 0
        ARG_filter_list = []
        for ii in range(len(ARG_list)):
            if ARG_list[ii][1] >= start_initial:
                ARG_filter_list.append((ARG_list[ii][0], ARG_list[ii][4], ARG_list[ii][5], ARG_list[ii][1], ARG_list[ii][2], ARG_list[ii][6], ARG_list[ii][7]))
                start_initial = ARG_list[ii][2]
                ii_keep = ii
            elif ARG_list[ii][1] < start_initial and ARG_list[ii][2] > start_initial and float(ARG_list[ii][2]-start_initial)/float(ARG_list[ii][2]-ARG_list[ii][1])>0.8 and float(start_initial- ARG_list[ii][1])/float(ARG_list[ii][0].split(":")[-1]) < 0.2:
                ARG_filter_list.append((ARG_list[ii][0], ARG_list[ii][4], ARG_list[ii][5], ARG_list[ii][1], ARG_list[ii][2], ARG_list[ii][6], ARG_list[ii][7]))
                start_initial = ARG_list[ii][2]
                ii_keep = ii
            else:
                if ARG_list[ii][3] > ARG_list[ii_keep][3]:
                    sss = (ARG_list[ii_keep][0],ARG_list[ii_keep][4], ARG_list[ii_keep][5], ARG_list[ii][1], ARG_list[ii][2], ARG_list[ii][6], ARG_list[ii][7])
                    if sss in ARG_filter_list:
                        ARG_filter_list.remove(sss)
                    ARG_filter_list.append((ARG_list[ii][0], ARG_list[ii][4], ARG_list[ii][5], ARG_list[ii][1], ARG_list[ii][2], ARG_list[ii][6], ARG_list[ii][7]))
                    start_initial = ARG_list[ii][2]
                    ii_keep = ii
                else:
                    continue
        ARG_filter_list_use = list(set(ARG_filter_list))
        for j in ARG_filter_list_use:
            if j not in ARG_list_get:
                ARG_list_get.append(j[0])
                ARG_list_get = list(set(ARG_list_get))
        ARG_location_filter_dict[key_data] = ARG_filter_list_use
    ARG_location_use_dict = {}
    for kk in ARG_location_filter_dict.items():
        ID_filter = kk[0]
        ARG_list = kk[1]
        Isolates_use = ID_filter[0]
        for i in ARG_list:
            ARG_location_use_dict.setdefault(Isolates_use,[]).append(i)
    ARG_list_get.sort()
    return ARG_location_use_dict, ARG_list_get

def normal_display():
    ARG_abstract_dict = {}
    ARG_location_dict = blast_filter()
    ARG_location_use_dict, ARG_list_get = ARG_filter(ARG_location_dict)
    out_file = open("data_results.out", "w")
    out_file.write("ID" + "\t")
    for line in ARG_list_get:
        out_file.write(line + "\t")
    out_file.write("\n")
    for item in ARG_location_use_dict.items():
        isolates_info = item[0]
        out_file.write(isolates_info + "\t")
        ARG_list = item[1]
        ARG_write_list = []
        ARG_write_dict = {}
        for ARG_have in ARG_list:
            ARG_write_list.append(ARG_have[0])
            ARG_write_dict.setdefault(ARG_have[0],[]).append((ARG_have[1], ARG_have[2], ARG_have[3], ARG_have[4], ARG_have[5], ARG_have[6]))
        for n in ARG_list_get:
            ARG_write = "0"
            if n in ARG_write_list:
                ARG_write_use = ARG_write_dict[n]
                for i in ARG_write_use:
                    if float(i[0]) >= 80.00 and float(i[1]) >= 80.00: # cutoff values for deretmining the mapped genes.
                        ARG_write = "1"
                        ARG_abstract_dict.setdefault(isolates_info,[]).append((n, i[4], i[2], i[3], i[5]))
            else:
                ARG_write = "0"
            out_file.write(str(ARG_write) + "\t")
        out_file.write("\n")
    out_file.close()
    return ARG_abstract_dict

def blast_seq_abstract(ARG_abstract_dict):
    rootdir1 = "./"
    for parent,dirnames,filenames in os.walk(rootdir1):
        for filename in filenames:
            aa = filename.split(".")[-1]
            if aa == "fna": # the extension of your genome file should be ".fna"
                bb = filename.split(".")[0]
                if bb in ARG_abstract_dict:
                    ARG_list = ARG_abstract_dict[bb]
                    file_data = open("./"+filename, "r")
                    out_file = open("./"+bb+ "_sc19_blast.fasta", "w") # the output abstract fasta file.
                    records_use = []
                    for record in SeqIO.parse(file_data,"fasta"):
                        ID = record.id
                        contig_len = len(record.seq)
                        for i in ARG_list:
                            if ID == i[1]:
                                ARG_ID = i[0]
                                F_R = i[4]
                                location_S = int(i[2])
                                location_E = int(i[3])
                                if F_R == "F": # F means the genes is Forward
                                    gene_seq = record.seq[location_S-1 :location_E]
                                    element_ID = bb + ":" + ARG_ID + ":" + F_R
                                    element_record = SeqRecord(gene_seq,id = element_ID, description = "")
                                    records_use.append(element_record)
                                elif F_R == "R": #R means the genes is Reverse
                                    gene_seq_ori = record.seq[location_S-1 :location_E]
                                    gene_seq = gene_seq_ori.reverse_complement()
                                    element_ID = bb + ":" + ARG_ID + ":" + F_R
                                    element_record = SeqRecord(gene_seq,id = element_ID, description = "")
                                    records_use.append(element_record)
                    SeqIO.write(records_use, out_file, "fasta")
                    out_file.close()
                    file_data.close()
                    
    return None


if __name__ == "__main__":
#    do_gene_blast()
#    ST_cal_display()
    ARG_abstract_dict = normal_display()
    blast_seq_abstract(ARG_abstract_dict)












