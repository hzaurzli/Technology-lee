# If you have got the core_gene_list_file.txt (a list of all core-genome genes, one gene per line), and all the mapped-genes-abstract files created by"Abstract_mapped_genes.py", just "cat" all the mapped-genes-abstract files into a single file named "all_abstract_genes.fna". and put the script together with above two files, run the script, you will got a file for each core gene that including the gene of all genomes.
# There is a little bug here. Sometimes there might be copies for a core gene in a genome, that will cause all the copies will be abstracted to the core genome file, and will cause in each core gene file the numbers of genomes are different. If you find this flaw, don't panic，this flaw will be polished in next steps, just keep running the script.




import re
import sys, getopt
import operator
from Bio import SeqIO



blast_info = open("./all_abstract_genes.fna", "r")
gene_info = open("./core_gene_list.txt", "r")


gene_list = []
for i in gene_info:
    gene_list.append(i.strip())


for record in SeqIO.parse(blast_info,"fasta"): #>083_2A:ADNHCAIG_00326:198:F
    ID_info = record.id.strip().split(":")
    Genome_ID = ID_info[0]
    Gene_ID = ID_info[1] + ":" + ID_info[2]
    if Gene_ID in gene_list:
        Gene_ID_use = Gene_ID.split(":")[0] + "-" + Gene_ID.split(":")[1]
        out_file = open("./" + Gene_ID_use + ".fa", "a")
        SeqIO.write(record, out_file,"fasta")
        out_file.close()
blast_info.close()

