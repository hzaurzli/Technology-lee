#-*-coding:utf-8-*-
import sys
import re
import sys, getopt
import operator
import os


## 读取 blast format 6 的内容
def blast_filter():
    ARG_location_dict = {}
    rootdir1 = "./"
    ARG_list = []
    for parent,dirnames,filenames in os.walk(rootdir1):
        for filename in filenames:
            cc = filename
            dd = cc.split(".")[0]
            ff = cc.split(".")[-1]
            if ff == "out":
                blast_info = open(filename, "r")
                for line in blast_info:
                    line_info = line.strip().split("\t")
                    Contig_ID_info = line_info[0]
                    ARG_ID = line_info[1]
                    ARG_len = line_info[1].split(":")[-1]
                    identical_percent = float(line_info[2])
                    align_length = int(line_info[3])
                    contig_start = int(line_info[6])
                    contig_end = int(line_info[7])
                    ARG_start = int(line_info[8])
                    ARG_end = int(line_info[9])
                    Score = float(line_info[11])
                    align_percent = float(align_length)/float(ARG_len)*100
                    key_use = (dd, Contig_ID_info)
                    ARG_location_dict.setdefault(key_use,[]).append((ARG_ID, contig_start, contig_end, Score, identical_percent, ARG_start, ARG_end))
                    ARG_list.append(ARG_ID)
                blast_info.close()
    ARG_list_filter = list(set(ARG_list))
    ARG_list_use = sorted(ARG_list_filter)
    return ARG_location_dict, ARG_list_use

## 挑选每一个sample 的 ref 被 align 到的区间,并从小打大排序
def calculate_coverage(ARG_location_dict):
    ARG_loctaion_cal_dict = {}
    for item in ARG_location_dict.items():
        ARG_key = item[0]
        strain_name = ARG_key[0]
        ARG_list = item[1]
        for i in ARG_list:
            ARG_ID = i[0]
            ARG_start = min(i[5], i[6])
            ARG_end = max(i[5], i[6])
            percentage = i[4]
            if percentage >= 80.00:
                key = (strain_name, ARG_ID)
                add_use = (ARG_start, ARG_end)
                ARG_loctaion_cal_dict.setdefault(key,[]).append(add_use)
    coverage_cal_dict = {}
    for ii in ARG_loctaion_cal_dict.items():
        strain_info = ii[0]
        mapped_list = ii[1]
        strain_ID = strain_info[0]
        ARG_len = strain_info[1].split(":")[-1]
        ARG_ID = strain_info[1]
        mapped_list.sort(key = operator.itemgetter(0,1))
        if len(mapped_list) == 1:
            Subject_cover_ca = float(mapped_list[0][1] - mapped_list[0][0]+1)
            # 覆盖的区间去重加和后除以 ref 某序列的长度,计算百分比
            Percent_Subject_ca = Subject_cover_ca/float(ARG_len)*100
        if len(mapped_list) > 1:
            Subject_cover_ca = int(mapped_list[0][1] - mapped_list[0][0]+1)
            num_ca = 1
            Big_contig_end_ca = int(mapped_list[0][1])
            while (num_ca <= len(mapped_list)-1):
                if int(mapped_list[num_ca][1]) <= int(Big_contig_end_ca):
                    num_ca += 1
                else:
                    if int(mapped_list[num_ca][0]) >= int(Big_contig_end_ca):
                        Subject_cover_ca = Subject_cover_ca +  int(mapped_list[num_ca][1]) - int(mapped_list[num_ca][0])+1
                        Big_contig_end_ca = int(mapped_list[num_ca][1])
                        num_ca += 1
                    else:
                        Subject_cover_ca = Subject_cover_ca + int(mapped_list[num_ca][1]) - int(Big_contig_end_ca)
                        Big_contig_end_ca = int(mapped_list[num_ca][1])
                        num_ca += 1
            # 覆盖的区间去重加和后除以 ref 某序列的长度,计算百分比
            Percent_Subject_ca = (float(Subject_cover_ca)/float(ARG_len))*100
        coverage_cal_dict.setdefault(strain_ID,[]).append((ARG_ID,Percent_Subject_ca))
    return coverage_cal_dict



def ARG_filter(ARG_location_dict):
#    list_mark = Mark_list()
    ARG_location_filter_dict = {}
    for item in ARG_location_dict.items():
        key_data = item[0]
        ARG_list = item[1]
        ARG_list.sort(key = operator.itemgetter(1))
        start_initial = 0
        ARG_filter_list = []
        for ii in range(len(ARG_list)):
#            if ARG_list[ii][0] in list_mark:
            if ARG_list[ii][1] >= start_initial:
                ARG_filter_list.append(ARG_list[ii][0])
                start_initial = ARG_list[ii][2]
            else:
                if ARG_list[ii][3] > ARG_list[ii - 1][3]:
                    sss = ARG_list[ii - 1]
                    if sss in ARG_filter_list:
                        ARG_filter_list.remove(sss)
                    ARG_filter_list.append(ARG_list[ii][0])
                    start_initial = ARG_list[ii][2]
                else:
                    continue

#            else:
#                ARG_filter_list.append(ARG_list[ii][0])
        ARG_filter_list_use = list(set(ARG_filter_list))
        ARG_location_filter_dict[key_data] = ARG_filter_list_use
    ARG_location_use_dict = {}
    for  kk in ARG_location_filter_dict.items():
        ID_filter = kk[0]
        ARG_list = kk[1]
        Isolates_use = ID_filter[0]
        for i in ARG_list:
            ARG_location_use_dict.setdefault(Isolates_use,[]).append(i)
    return ARG_location_use_dict



def normal_display():
    ARG_location_dict, ARG_list_use = blast_filter()
    coverage_cal_dict = calculate_coverage(ARG_location_dict)
    out_file = open("data_results.out", "w")
    out_file_2 = open("best_match.out", "w")
    out_file.write("ID" + "\t")
    for line in ARG_list_use:
        out_file.write(line + "\t")
    out_file.write("\n")
    for item in coverage_cal_dict.items():
        isolates_info = item[0]
        out_file.write(isolates_info + "\t")
        ARG_list = item[1]
        ARG_list.sort(key = operator.itemgetter(1), reverse=True)
        out_file_2.write(isolates_info + "\t" + str(ARG_list[0][0])+ "\t" + str(ARG_list[0][1]) + "\n")
        ARG_write_list = []
        ARG_write_dict = {}
        ARG_write = "0"
        for ll in ARG_list:
            ARG_have = ll[0]
            ARG_write_list.append(ARG_have)
            ARG_write_dict[ARG_have] = ll[1]
        for n in ARG_list_use:
            if n in ARG_write_list and ARG_write_dict[n] >= 60:
                ARG_write = ARG_write_dict[n]
            else:
                ARG_write = "0"
            out_file.write(str(ARG_write) + "\t")
        out_file.write("\n")
    out_file.close()
    out_file_2.close()
    return None



if __name__ == "__main__":
    normal_display()













