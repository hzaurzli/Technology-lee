# This is the script for abstract mapped genes from blast files.(the extension is .blast. You can change it in the script...do not forget the blast output formate is 6).
# Just put all your genome files (the extension is .fna) under a folder named "fna", and Reference genome database in a folder named "Ref_database", and put the two folders as well as the script under a same folder. Run the script, you will got the algnment file and SNP matrix at the same time, but just make sure you have right install the alignment program "clustalo" and snp abstract program "snp-dists", and you have to change the associated paths of the programs in the script.

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

def do_gene_blast(): #
    rootdir1 = "./fna"
    rootdir2 = "./Ref_database/"
    for parent2,dirnames2,filenames2 in os.walk(rootdir2):
        for filename2 in filenames2:
            qq = filename2
            zz = qq.split(".")[-1]
            if zz == "nhr":
                bb = qq.split(".")[0]
    for parent,dirnames,filenames in os.walk(rootdir1):
        for filename in filenames:
            cc = filename
            dd = cc.split(".")[0]
            ff = cc.split(".")[-1]
            if ff == "fna":
                os.system("/home/rzli/software/ncbi-blast-2.13.0+/bin/blastn -db " +  rootdir2 + bb + " -query ./fna/" + cc + " -out ./fna/" + dd + ".blast -outfmt 6 -evalue 1e-5 -max_target_seqs 100000")
    return None

def blast_filter():
    ARG_location_dict = {}
    rootdir1 = "./fna"
    for parent,dirnames,filenames in os.walk(rootdir1):
        for filename in filenames:
            cc = filename
            dd = cc.split(".")[0]
            ff = cc.split(".")[-1]
            if ff == "blast" : # the extension of the blast files
                blast_info = open("./fna/" + filename, "r")
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
                    if float(align_percent) >= 80.00 and float(identical_percent) >= 80.00:
                        ARG_location_dict.setdefault(key_use,[]).append((ARG_ID, contig_start, contig_end, Score, identical_percent, align_percent, Contig_ID_info, F_R))
                blast_info.close()
    return ARG_location_dict

def ARG_filter(ARG_location_dict):
    ARG_location_filter_dict = {}
    for item in ARG_location_dict.items():
        key_data = item[0]
        ARG_list = item[1]
        ARG_list.sort(key = operator.itemgetter(1))
        start_initial = 0
        ii_keep = 0
        ARG_filter_list = []
        for ii in range(len(ARG_list)):
            if ARG_list[ii][1] >= start_initial:
                ARG_filter_list.append((ARG_list[ii][0], ARG_list[ii][4], ARG_list[ii][5], ARG_list[ii][1], ARG_list[ii][2], ARG_list[ii][6], ARG_list[ii][7], ARG_list[ii][3]))
                start_initial = ARG_list[ii][2]
                ii_keep = ii
            elif ARG_list[ii][1] < start_initial and ARG_list[ii][2] > start_initial and float(ARG_list[ii][2]-start_initial)/float(ARG_list[ii][2]-ARG_list[ii][1])>0.8 and float(start_initial- ARG_list[ii][1])/float(ARG_list[ii][0].split(":")[-1]) < 0.2:
                ARG_filter_list.append((ARG_list[ii][0], ARG_list[ii][4], ARG_list[ii][5], ARG_list[ii][1], ARG_list[ii][2], ARG_list[ii][6], ARG_list[ii][7],ARG_list[ii][3]))
                start_initial = ARG_list[ii][2]
                ii_keep = ii
            else:
                if ARG_list[ii][3] > ARG_list[ii_keep][3]:
                    sss = (ARG_list[ii_keep][0],ARG_list[ii_keep][4], ARG_list[ii_keep][5], ARG_list[ii][1], ARG_list[ii][2], ARG_list[ii][6], ARG_list[ii][7],ARG_list[ii][3])
                    if sss in ARG_filter_list:
                        ARG_filter_list.remove(sss)
                    ARG_filter_list.append((ARG_list[ii][0], ARG_list[ii][4], ARG_list[ii][5], ARG_list[ii][1], ARG_list[ii][2], ARG_list[ii][6], ARG_list[ii][7], ARG_list[ii][3]))
                    start_initial = ARG_list[ii][2]
                    ii_keep = ii
                else:
                    continue
        ARG_filter_list_use = list(set(ARG_filter_list))
        ARG_location_filter_dict[key_data] = ARG_filter_list_use
    ARG_location_use_dict = {}
    for kk in ARG_location_filter_dict.items():
        ID_filter = kk[0]
        ARG_list = kk[1]
        Isolates_use = ID_filter[0]
        for i in ARG_list:
            ARG_location_use_dict.setdefault(Isolates_use,[]).append(i)
    return ARG_location_use_dict

def Core_genome_cal(ARG_location_use_dict):
    Gene_dict = {}
    Genome_gene_have_dict = {}
    isolates_number = len(ARG_location_use_dict)
    Gene_number_list = []
    core_gene_list = []
    for item in ARG_location_use_dict.items():
        isolate_ID = item[0]
        Gene_list = item[1]
        for i in Gene_list:
            Gene_ID = i[0]
            Gene_identical = i[1]
            Gene_aln_percent = i[2]
            Gene_contig_loc_star = i[3]
            Gene_contig_loc_end = i[4]
            Contig_info = i[5]
            Gene_F_R = i[6]
            Gene_Score = i[7]
            Gene_dict.setdefault(Gene_ID, []).append(isolate_ID)
            Genome_gene_have_dict.setdefault(isolate_ID, []).append(Gene_ID)
    for j in Genome_gene_have_dict.items():
        isolate_name = j[0]
        Gene_list = j[1]
        Gene_list_filter = list(set(Gene_list))
        Gene_number_list.append((isolate_name, len(Gene_list_filter)))
    Gene_number_list.sort(key = operator.itemgetter(1))
    out_file = open("./Isolates_gene_number.out", "w")
    out_file2 = open("./core_gene_list.out", "w")
    out_file.write("Genome_ID" + "\t" + "Gene_number" + "\n")
    print("The genome with least gene number is: ", Gene_number_list[0][0], ":", Gene_number_list[0][1])
    for k in Gene_number_list:
        out_file.write(str(k[0]) + "\t" + str(k[1]) + "\n")
    out_file.close()
    for i in Gene_dict.items():
        Gene_name = i[0]
        Genome_list = i[1]
        Genome_list_filter = list(set(Genome_list))
        if len(Genome_list_filter) == isolates_number:
            core_gene_list.append(Gene_name)
    for m in core_gene_list:
        out_file2.write(m + "\n")
    print("The number of core genes is : ", len(core_gene_list))
    print("Do you want to remove some genomes? ")
    out_file2.close()
    out_file.close()
    return core_gene_list
        
def Gene_abstract_determine(ARG_location_use_dict, core_gene_list):
    ARG_abstract_dict = {}
    ARG_location_dict = blast_filter()
    for item in ARG_location_use_dict.items():
        Gene_filter_dict = {}
        Gene_list_filter = []
        isolates_info = item[0]
        Gene_list = item[1]
        for i in Gene_list:
            Gene_ID = i[0]
            Gene_identical = i[1]
            Gene_aln_percent = i[2]
            Gene_contig_loc_star = i[3]
            Gene_contig_loc_end = i[4]
            Contig_info = i[5]
            Gene_F_R = i[6]
            Gene_Score = i[7]
            if Gene_ID in core_gene_list:
                Gene_filter_dict.setdefault(Gene_ID, []).append(i)
        for item_i in Gene_filter_dict.items():
            Gene_use = item_i[0]
            Gene_count_list = item_i[1]
            if len(Gene_count_list) == 1:
                Gene_list_filter.append(Gene_count_list[0])
            else:
                Gene_count_list.sort(key = operator.itemgetter(-1))
                Gene_use_final = Gene_count_list[-1]
                Gene_list_filter.append(Gene_use_final)
        ARG_abstract_dict[isolates_info] = Gene_list_filter
    return ARG_abstract_dict

def blast_seq_abstract(ARG_abstract_dict):
    rootdir1 = "./fna"
    blast_info = open("./All_isolates_core_genes.fna", "w")
    records_use = []
    for parent,dirnames,filenames in os.walk(rootdir1):
        for filename in filenames:
            aa = filename.split(".")[-1]
            if aa == "fna": # the extension of your genome file should be ".fna"
                bb = filename.split(".")[0]
                if bb in ARG_abstract_dict:
                    ARG_list = ARG_abstract_dict[bb]
                    file_data = open("./fna/"+filename, "r")
                    for record in SeqIO.parse(file_data,"fasta"):
                        ID = record.id
                        contig_len = len(record.seq)
                        for i in ARG_list:
                            if ID == i[5]:
                                ARG_ID = i[0]
                                F_R = i[6]
                                location_S = int(i[3])
                                location_E = int(i[4])
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
                    file_data.close()
    SeqIO.write(records_use,  blast_info, "fasta")
    blast_info.close()
    return records_use


def creat_aln_file(records_use, core_gene_list):
    gene_list = []
    core_gene_list_use = sorted(core_gene_list)
    for record in records_use: #>083_2A:ADNHCAIG_00326:198:F
        ID_info = record.id.strip().split(":")
        Genome_ID = ID_info[0]
        Gene_ID = ID_info[1] + ":" + ID_info[2]
        if Gene_ID in core_gene_list_use:
            Gene_ID_use = Gene_ID.split(":")[0] + "-" + Gene_ID.split(":")[1]
            out_file = open("./" + Gene_ID_use + ".corefa", "a")
            SeqIO.write(record, out_file,"fasta")
            out_file.close()
    return None
    


def do_alignment():
    rootdir1 = "./"
    for parent,dirnames,filenames in os.walk(rootdir1):
        for filename in filenames:
            aa = filename.split(".")[-1]
            if aa == "corefa":
                os.system("/home/rzli/miniconda3/bin/clustalo -i " + filename + " -o " + filename + ".aln -v")
    return None

def connect_aln_file():
    os.system("cat ./*.aln > All_core_genes_aln.fasta")
    blast_info = open("./All_core_genes_aln.fasta", "r")
    out_file = open("./All_SS_core_genes.aln", "w")
    genome_seq_dict = {}
    for record in SeqIO.parse(blast_info,"fasta"):
        ID = record.id
        isolate_use = ID.strip().split(":")[0]
        core_gene = ID.strip().split(":")[1]
        Seq = record.seq
        Len_use = len(Seq)
        Len_key = (Len_use, core_gene)
        genome_seq_dict.setdefault(isolate_use, []).append((core_gene, Seq))
    record_list = []
    for item in genome_seq_dict.items():
        Seq_use = ""
        genome_ID = item[0]
        Seq_list = item[1]
        Seq_list_use = sorted(Seq_list)
        for jj in Seq_list_use:
            Seq_use = Seq_use + jj[1]
        Seq_record = SeqRecord(Seq_use, id = genome_ID, description = "")
        record_list.append(Seq_record)
    SeqIO.write(record_list, out_file,"fasta")
    out_file.close()
    blast_info.close()
    return None

if __name__ == "__main__":
    do_gene_blast()
    ARG_location_dict = blast_filter()
    ARG_location_use_dict = ARG_filter(ARG_location_dict)
    core_gene_list = Core_genome_cal(ARG_location_use_dict)
    ARG_abstract_dict = Gene_abstract_determine(ARG_location_use_dict, core_gene_list)
    records_use = blast_seq_abstract(ARG_abstract_dict)
    creat_aln_file(records_use, core_gene_list)
    do_alignment()
    connect_aln_file()
    os.system("/home/rzli/miniconda3/envs/snp-dists/bin/snp-dists All_SS_core_genes.aln > All_SS_SNP.tsv")












