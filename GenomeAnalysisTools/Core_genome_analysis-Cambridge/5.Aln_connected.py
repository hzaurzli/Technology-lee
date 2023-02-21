# If you have got all the alignment files of core genes. Just "cat" all the files into a single file named "All_core_genes_aln.fasta", and put the script with this file under a same folder and running the script, you will got the final core genome alignment file of all your genomes.
# I have mentioned a flaw in the previous script, it was that in the alignment file of each core gene, there might be copies of the core gene sequences from a same genome. I have to choose the best one and remove the other copies. I just count the gap "-" in the genes, the one with least "-" will be kept (I know the number of mismatch gene should also be included, but in the further, I will fix it in the first step to just compare the blast score to choose the best one.)


import re
import os
import sys, getopt
import operator
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

blast_info = open("./All_core_genes_aln.fasta", "r")
out_file = open("./All_SS_core_genes.aln", "w")
#out_check = open("./Length_check.out", "w")
genome_seq_dict = {}
#length_check_genes_dict = {}
for record in SeqIO.parse(blast_info,"fasta"):
    ID = record.id
    isolate_use = ID.strip().split(":")[0]
    core_gene = ID.strip().split(":")[1]
    Seq = record.seq
    Len_use = len(Seq)
    Len_key = (Len_use, core_gene)
#    length_check_genes_dict.setdefault(Len_key, []).append(isolate_use)
    genome_seq_dict.setdefault(isolate_use, []).append((core_gene, Seq))

#for kkk in length_check_genes_dict.items():
#   kkk_key = kkk[0]
#   kkk_list = kkk[1]
#   out_check.write(kkk_key[1] + "\t" + str(kkk_key[0]))
#   for zzz in kkk_list:
#       out_check.write("\t" +zzz)
#   out_check.write("\n")

genome_seq_filter_dict = {}
for fff in genome_seq_dict.items():
    fff_ID = fff[0]
    fff_Seq_list = fff[1]
    filter_dict = {}
    for ff in fff_Seq_list:
        filter_dict.setdefault(ff[0], []).append(ff[1])
    for ff_item in filter_dict.items():
        ff_core_gene_ID = ff_item[0]
        ff_core_gene_list = ff_item[1]
        if len(ff_item[1]) == 1:
            genome_seq_filter_dict.setdefault(fff_ID, []).append((ff_core_gene_ID, ff_core_gene_list[0]))
        else:
            n_count_list = []
            for ww in ff_core_gene_list:
                n_count = ww.count("-")
                n_count_list.append((n_count, ww))
            n_count_list_use = sorted(n_count_list)
            genome_seq_filter_dict.setdefault(fff_ID, []).append((ff_core_gene_ID, n_count_list_use[0][1]))
    
record_list = []
for item in genome_seq_filter_dict.items():
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
out_check.close()


