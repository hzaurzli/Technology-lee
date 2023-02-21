# This is the script for Streptococcus suis serotyping, assembled contig file input (nice!). Yes, it can distiguish serotype 1 from serotype 14 and serotype 2 from serotype 1/2. The cutoff value for the identification of a serotyping-determing gene is 90% identity and 90% coverage, you can change it in the script.

#How to use the script? Just put all the serotyping blast result files (the extension is .Sero, if you dislike it you can change it in the script...do not forget the blast output formate is 6) of your isolates with the script under the same folder, and run the script. You will got the result file named SS_serotype_results.out.

# Maybe some weird results will in your results, for example, something like "cpsK:Serotype_10:Serotype_1", some serotype linked by ":", just means that there are more than one specific serotyping-determing genes were identified in your genome, maybe your genome data is heterogeneity.

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
from Bio import pairwise2
from Bio.Seq import translate

def blast_filter():
    ARG_location_dict = {}
    isolates_list = []
    rootdir1 = "./"
    for parent,dirnames,filenames in os.walk(rootdir1):
        for filename in filenames:
            cc = filename
            dd = cc.split(".")[0]
            ff = cc.split(".")[-1]
            if ff == "Sero": # the extension of your blast results file
                blast_info = open(filename, "r")
                isolates_list.append(dd)
                for line in blast_info:
                    line_info = line.strip().split("\t")
                    Contig_ID_info = line_info[0]
                    ARG_ID = line_info[1].split(":")[0]
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
                    key_use = (dd)
                    if float(align_percent) > 50.00:
                        ARG_location_dict.setdefault(key_use,[]).append((ARG_ID, contig_start, contig_end, Score, identical_percent, align_percent, Contig_ID_info, F_R))
                blast_info.close()
    return ARG_location_dict, isolates_list

def ARG_filter(ARG_location_dict, isolates_list):
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
        for i in ARG_list:
            if float(i[1]) >= 90.00 and float(i[2]) >= 90.00: # the cutoff values for the identification of a serotyping-determing gene
                ARG_location_use_dict.setdefault(ID_filter,[]).append(i)
    for jj in isolates_list:
        if jj not in ARG_location_use_dict:
            ARG_location_use_dict[jj] = "Nontypeable"
    return ARG_location_use_dict, ARG_list_get


def further_typing(gene_protein, serotype_initial):
    final_serotype = "Waring"
    ref_seq = Seq("MINISIIVPIYNVEQYLSKCINSIVNQTYKHIEILLVNDGSTDNSEEICLAYAKKDSRIRYFKKENGGLSDARNYGISRAKGDYLAFIDSDDFIHSEFIQRLHEAIERENALVAVAGYDRVDASGHFLTAEPLPTNQAVLSGRNVCKKLLEADGHRFVVAWNKLYKKELFEDFRFEKGKIHEDEYFTYRLLYELEKVAIVKECLYYYVDRENSITTSSMTDHRFHCLLEFQNERMDFYESRGDKELLLECYRSFLAFAVLFLGKYNHWLSKQQKKLLQTLFRIVYKQLKQNKRLALLMNAYYLVGCLHLNFSVFLKTGKDKIQERLRRSESSTR")
    Query_seq = Seq(gene_protein)
    Alignments_use = pairwise2.align.localxx(ref_seq, Query_seq)
    Seq_A = Alignments_use[-1].seqA
    Seq_B = Alignments_use[-1].seqB
    n = 0
    m = 0
    for i in Seq_A:
        if i != "-":
            m = m + 1
            n = n + 1
        else:
            n = n + 1
        if m == 161:
            break
    Check_A = Seq_A[n-1]
    Check_B = Seq_B[n-1]
    print(Check_A, Check_B, serotype_initial)
    n = 0
    m = 0
    if Check_A != "W":
        print("Wrong!")
    if Check_B == "W" and serotype_initial == "Serotype_1":
        final_serotype = "Serotype_14"
    if Check_B == "C" and serotype_initial == "Serotype_1":
        final_serotype = "Serotype_1"
    if Check_B == "W" and serotype_initial == "Serotype_2":
        final_serotype = "Serotype_2"
    if Check_B == "C" and serotype_initial == "Serotype_2":
        final_serotype = "Serotype_1/2"
    if Check_B != "C" and Check_B != "W":
        around_seq = Seq_B[n-5]+ Seq_B[n-4]+Seq_B[n-3]+Seq_B[n-2]+Seq_B[n-1]+Seq_B[n]+Seq_B[n+1]+Seq_B[n+2]+Seq_B[n+3]
        final_serotype = around_seq
    return final_serotype

def sero_typing():
    ARG_location_dict, isolates_list = blast_filter()
    ARG_location_use_dict, ARG_list_get = ARG_filter(ARG_location_dict, isolates_list)
    out_file = open("SS_serotype_results.out", "w")
    out_file.write("Genome_ID" + "\t" + "Serotype" + "\n")
    for i in isolates_list:
        Serotype_results_info = ARG_location_use_dict[i]
        if Serotype_results_info == "Nontypeable":
            out_file.write(i + "\t" + "Nontypeable" + "\n")
        elif len(Serotype_results_info) == 1 and Serotype_results_info[0][0] != "Serotype_1" and Serotype_results_info[0][0] != "Serotype_2" and Serotype_results_info[0][0] != "cpsK":
            out_file.write(i + "\t" + Serotype_results_info[0][0] + "\n")
        elif len(Serotype_results_info) == 1 and Serotype_results_info[0][0] == "cpsK":
            out_file.write(i + "\t" + "Nontypeable" + "\n")
        elif len(Serotype_results_info) == 1 and (Serotype_results_info[0][0] == "Serotype_1" or Serotype_results_info[0][0] == "Serotype_2"):
            if Serotype_results_info[0][0] == "Serotype_1":
                out_file.write(i + "\t" + "Serotype_1&14"  + "\n")
            elif Serotype_results_info[0][0] == "Serotype_2":
                out_file.write(i + "\t" + "Serotype_2&1/2"  + "\n")
        elif len(Serotype_results_info) == 2:
            zz_list = []
            for zz in Serotype_results_info:
                zz_list.append(zz[0])
            if "cpsK" in zz_list and ("Serotype_1" not in zz_list) and ("Serotype_2" not in zz_list):
                    out_file.write(i + "\t" + zz[0] + "\n")
            elif "cpsK" in zz_list and (("Serotype_1" in zz_list) or ("Serotype_2" in zz_list)) :
                final_serotype = "Weird!"
                rootdir1 = "./" #将所有基因组fasta(后缀名需为.fna)序列放入当前文件夹
                for parent,dirnames,filenames in os.walk(rootdir1):
                    for filename in filenames:
                        aa = filename.split(".")[-1]
                        if aa == "fna":
                            bb = filename.split(".")[0]
                            if bb == i:
                                for kk in Serotype_results_info:
                                    if kk[0] == "cpsK":
                                        file_data = open("./"+filename, "r")
                                        records_use = []
                                        for record in SeqIO.parse(file_data,"fasta"):
                                            ID = record.id
                                            if ID == kk[5]:
                                                contig_len = len(record.seq)
                                                ARG_ID = kk[0]
                                                F_R = kk[6]
                                                location_S = int(kk[3])
                                                location_E = int(kk[4])
                                                if F_R == "F":
                                                    gene_seq = record.seq[location_S-1 :location_E]
                                                    frame_list = []
                                                    for frame_start in range(3):
                                                        frame = translate(gene_seq[frame_start:],table=11)
                                                        Non_count = frame.count("*")
                                                        frame_list.append((Non_count,frame))
                                                    frame_list_use = sorted(frame_list)
                                                    gene_protein = frame_list_use[0][1]
                                                    if "Serotype_1" in zz_list:
                                                        final_serotype = further_typing(gene_protein, "Serotype_1")
                                                    elif "Serotype_2" in zz_list:
                                                        final_serotype = further_typing(gene_protein, "Serotype_2")
                                                    else:
                                                        final_serotype = "Waring!"
                                                elif F_R == "R":
                                                    gene_seq_ori = record.seq[location_S-1 :location_E]
                                                    gene_seq = gene_seq_ori.reverse_complement()
                                                    for frame_start in range(3):
                                                        frame = translate(gene_seq[frame_start:],table=11)
                                                        Non_count = frame.count("*")
                                                        frame_list.append((Non_count,frame))
                                                    frame_list_use = sorted(frame_list)
                                                    gene_protein = frame_list_use[0][1]
                                                    if "Serotype_1" in zz_list:
                                                        final_serotype = further_typing(gene_protein, "Serotype_1")
                                                    elif "Serotype_2" in zz_list:
                                                        final_serotype = further_typing(gene_protein, "Serotype_2")
                                                    else:
                                                        final_serotype = "Waring!"
                                        file_data.close()
                out_file.write(i + "\t" + final_serotype + "\n")
            else:
                out_file.write(i + "\t" + zz[0]+":"+zz[1]+ "\n")
        else:
            Weird_list = []
            for kkk in Serotype_results_info:
                Weird_list.append(kkk[0])
            Write_Werid = ":".join(Weird_list)
            out_file.write(i + "\t" + Write_Werid + "\n")
    out_file.close()
    return None


if __name__ == "__main__":
#    do_gene_blast()
#    ST_cal_display()
    sero_typing()












