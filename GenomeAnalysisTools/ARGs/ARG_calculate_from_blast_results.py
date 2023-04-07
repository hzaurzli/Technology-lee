#This is the script for calculating ARG or VFs genes from blast results, and try to find new alleles. The blast output formate should be 6. The ARG/VF blast database should be modified by adding the length of genes to the end of the gene ID and saperated with ":"
# Put the blast results file under a folder named /ARG_blast/, and put the script with the /ARG_blast/ folder under a same folder. Run the script.
# This is just a home used script, not be polished, very personal style and no debugging, maybe not suitable for public using. Writtern by Geng ZOU,on 2022-11-02, if you have any questions, contact me. zougeng19900918@126.com
#-*-coding:utf-8-*-
import sys
import re
import sys, getopt
import operator
import os

def blast_filter():
    ARG_location_dict = {}
    for parent,dirnames,filenames in os.walk(rootdir1):
        for filename in filenames:
            cc = filename
            dd = cc.split(".")[0]
            ff = cc.split(".")[-1]
            if ff == "ARG" :
                blast_info = open(rootdir1 + filename, "r")
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
                    align_percent = '%.2f' %(float(align_length)/float(ARG_len)*100)
                    key_use = (dd, Contig_ID_info)
                    if float(align_percent) > 50.00: # this is the cutoff value for finding new alleles, the alignment length coverage larger than 50%
                        ARG_location_dict.setdefault(key_use,[]).append((ARG_ID, contig_start, contig_end, Score, identical_percent, align_percent))
                blast_info.close()
    print(ARG_location_dict)
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
                ARG_filter_list.append((ARG_list[ii][0], ARG_list[ii][4], ARG_list[ii][5]))
                start_initial = ARG_list[ii][2]
                ii_keep = ii
            elif ARG_list[ii][1] < start_initial and ARG_list[ii][2] > start_initial and float(ARG_list[ii][2]-start_initial)/float(ARG_list[ii][2]-ARG_list[ii][1])>0.8 and float(start_initial- ARG_list[ii][1])/float(ARG_list[ii][0].split(":")[-1]) < 0.2:
                ARG_filter_list.append((ARG_list[ii][0], ARG_list[ii][4], ARG_list[ii][5]))
                start_initial = ARG_list[ii][2]
                ii_keep = ii
            else:
                if ARG_list[ii][3] > ARG_list[ii_keep][3]:
                    sss = (ARG_list[ii_keep][0],ARG_list[ii_keep][4], ARG_list[ii_keep][5])
                    if sss in ARG_filter_list:
                        ARG_filter_list.remove(sss)
                    ARG_filter_list.append((ARG_list[ii][0], ARG_list[ii][4], ARG_list[ii][5]))
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
    ARG_location_dict = blast_filter()
    ARG_location_use_dict, ARG_list_get = ARG_filter(ARG_location_dict)
    out_file = open(rootdir1+"data_results.out", "w")
    out_file_2 = open(rootdir1+"new_alleles.out", "w")
#    out_file.write("ID" + "\t")
#    for line in ARG_list_get:
#        out_file.write(line + "\t")
#    out_file.write("\n")
    for item in ARG_location_use_dict.items():
        isolates_info = item[0]
        out_file.write(isolates_info + "\t")
        out_file_2.write(isolates_info + "\t")
        ARG_list = item[1]
        ARG_write_list = []
        ARG_write_dict = {}
        for ARG_have in ARG_list:
            ARG_write_list.append(ARG_have[0])
            ARG_write_dict.setdefault(ARG_have[0],[]).append((ARG_have[1], ARG_have[2]))
        for n in ARG_list_get:
            ARG_write = "0"
            m = 0
            if n in ARG_write_list:
                ARG_write_use = ARG_write_dict[n]
                for i in ARG_write_use:
                    if float(i[0]) >= 80.00 and float(i[1]) >= 80.00: # this is the cutoff value for determining an ARG/VF gene
                        m = m + 1
                        ARG_write = m
                        if m > 0 :
                            out_file.write(n + ":" + str(i[0]) + ":" + str(i[1]) + "\t")
                    else:
                        out_file_2.write(n + ":" + str(i[0]) + ":" + str(i[1]) + "\t")
            else:
                ARG_write = "0"
#            out_file.write(str(ARG_write) + "\t")
        out_file.write("\n")
        out_file_2.write("\n")
    out_file.close()
    return None

if __name__ == "__main__":
    rootdir1 = "/home/rzli/ssuis/LSSP/faa/" 
    normal_display()
