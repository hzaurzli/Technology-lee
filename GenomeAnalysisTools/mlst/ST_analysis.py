import re
import sys, getopt
from itertools import islice
import operator


opts, args = getopt.getopt(sys.argv[1:], "i:")
blast_info = ""
for op, value in opts:
    if op =="-i":
        blast_info = open(value, "r")
        file_name = str(value)
file_name_partten = re.compile(r"^(.*).filter")
n = file_name_partten.search(file_name)
Strain_index = n.group(1)
Strain_type = "1"
ST_gene_list = []
ST_info_dict = {}
Warning = ""
info_file = open("./SS_ST_type_20221007.txt", "r")
out_file = open("./SS_ST.result", "a")
for line_ST_old in info_file:
    line_ST = line_ST_old.replace("\xc2\xa0", "")
    ST_info = line_ST.strip().split("\t")
    ST_type = str(ST_info[0])
    ST_str = ":".join(ST_info[1:8])
    key = str(ST_str)
    ST_info_dict[key] = (ST_type)
for line in blast_info:
    info = line.strip().split("\t")
    ST_gene = str(info[1])
    ST_gene_list.append(ST_gene)
v = []
for i in ST_gene_list:
    q = filter(str.isalpha, i)
    v.append(q)
v_f = list(set(v))
if len(v_f) != len(v):
    Warning = "Warning"
aroA = 0
cpn60 = 0
dpr = 0
gki = 0
mutS = 0
recA = 0
thrA = 0
ST_result = "new"
for item in ST_gene_list:
    if item[0] == "a" and item[1] == "r" and item[2] == "o":
        item_info = item.split(":")[0]
        item_info_get = item_info.split("_")[1]
        aroA = item_info_get
    elif item[0] == "c" and item[1] == "p" and item[2] == "n":
        item_info = item.split(":")[0]
        item_info_get = item_info.split("_")[1]
        cpn60 = item_info_get
    elif item[0] == "d" and item[1] == "p":
        item_info = item.split(":")[0]
        item_info_get = item_info.split("_")[1]
        dpr = item_info_get
    elif item[0] == "g" and item[1] == "k":
        item_info = item.split(":")[0]
        item_info_get = item_info.split("_")[1]
        gki =  item_info_get
    elif item[0] == "m":
        item_info = item.split(":")[0]
        item_info_get = item_info.split("_")[1]
        mutS =  item_info_get
    elif item[0] == "r":
        item_info = item.split(":")[0]
        item_info_get = item_info.split("_")[1]
        recA = item_info_get
    elif item[0] == "t":
        item_info = item.split(":")[0]
        item_info_get = item_info.split("_")[1]
        thrA = item_info_get
ST_test = str(aroA) + ":" + str(cpn60) + ":" + str(dpr) + ":" + str(gki) + ":" + str(mutS) + ":" + str(recA) + ":" + str(thrA)
ST_use = str(ST_test)
if ST_info_dict.has_key(ST_use)==True:
    ST_result = ST_info_dict[ST_test]
out_file.write(file_name + "\t" + ST_result + "\t" + Warning + "\t" + ST_use + "\n")

out_file.close()
info_file.close()
blast_info.close()
                         


