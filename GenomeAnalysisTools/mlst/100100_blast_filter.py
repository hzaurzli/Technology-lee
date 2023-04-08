import re
import sys, getopt

opts, args = getopt.getopt(sys.argv[1:], "i:o:")
blast_info = ""
out_file = ""
for op, value in opts:
    if op =="-i":
        blast_info = open(value, "r")
    elif op =="-o":
        out_file = open(value, "w")

for line in blast_info:
    info = line.strip().split("\t")
    PerCentID = float(info[2])
    if PerCentID == 100.00:
        AlignLen = int(info[3])
        Subject = str(info[1])
        Subject_pattern = re.compile(r"^.*\:(.*)")
        m = Subject_pattern.search(Subject)
        SubjectLen = m.group(1)
        Subject_figure = float(SubjectLen)
        Percentlen = float(AlignLen/Subject_figure)*100
        if Percentlen == 100.00:
            out_file.write("\t".join([info[0],info[1],info[2],info[3],info[4],info[5],info[6],info[7],info[8],info[9],info[10],info[11]])+"\n")

blast_info.close()
out_file.close()


        

