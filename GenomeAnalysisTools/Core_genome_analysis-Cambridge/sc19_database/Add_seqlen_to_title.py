# This little script is for adding the sequence length to the ID behind a colon ":" in a fasta file
# To run the script: python Add_seqlen_to_title.py -i YOUR_FASTA_FILE.fasta -o YOUR_OUTPUT_FILE.fasta


import re
import sys, getopt
import operator
from Bio import SeqIO

opts, args = getopt.getopt(sys.argv[1:], "hi:o:")
blast_info = ""
out_file = ""
for op, value in opts:
    if op =="-o":
        out_file = open(value, "w")
    if op =="-i":
        blast_info = open(value, "r")
    elif op == "-h":
        usage()
        sys.exit()

for record in SeqIO.parse(blast_info,"fasta"):
    Len = str(len(record.seq))
    ID = record.id
    record.id = ID+":"+Len
    record.description = ''
    SeqIO.write(record, out_file,"fasta")
out_file.close
