import os
import sys
from collections import defaultdict
import argparse

def mergedic(file):
    file = open(file,'r')
    dict_1 = {}
    for line in file.readlines():
        line = line.replace("\n", "")
        k = line.split('\t')[0]
        v = int(line.split('\t')[1])
        dict_1[k] = v
    return dict_1

def cbinddic(dic_1):
    key_1 = [key for key in dic_1]
    dic = {}
    for key in key_1:
        dic.setdefault(key,[]).append(dic_1[key])
    return dic
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Multi-value dictionary")
    parser.add_argument("-a", "--a", required=True, type=str, help="First file")
    parser.add_argument("-b", "--b", required=True, type=str, help="Files dirction")
    Args = parser.parse_args()
    file_1 = Args.a   #"/home/rzli/fish/zreb/0000/SEX0108-1-1_FRRB19H001013-1a.htseq"
    a_dic = mergedic(file_1)
    b_dic = cbinddic(a_dic)
    for i in os.listdir(Args.b):
        fullpath = Args.b + i
        if fullpath == file_1:
            continue
        else:
            a_file = mergedic(Args.b + i)
            for key in a_file:
                b_dic.setdefault(key, []).append(a_file[key])
    print(b_dic)
