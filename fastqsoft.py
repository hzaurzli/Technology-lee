#from __future__ import division
import argparse
import os.path

parser = argparse.ArgumentParser()
parser.add_argument('-fa',action='store_true', dest='fasta',default=False, help='Check out the sequance')
parser.add_argument('-nt',action='store_true', dest='ntCounts',default=False, help='Calculate Nucleotide count')
parser.add_argument('-gc',action='store_true', dest='gcContent',default=False, help='Calculate gc content ')
parser.add_argument('-r', action='store_true', dest='rnaSeq', default=False, help='Calculate Transcription sequence')

parser.add_argument('-i',dest='filename',type=argparse.FileType('r'),help='Input file') #读取文件
args = parser.parse_args()


fasta = {}
# with open(filenames) as file:
for line in args.filename:
    if line.startswith(">"):
        name = line[1:].rstrip()
        fasta[name] = ''
        continue
    fasta[name] += line.rstrip().upper()
# print(fasta)


for name in fasta.keys():
    ntCounts = []
    seq = fasta[name]
    for nt in ['A', 'C', 'G', 'T']:
        ntCounts.append(seq.count(nt))


    total = len(seq)
    gcCount = seq.count('G') + seq.count('C')
    gcContent = format(float(gcCount / total * 100),'.6f')


    rnaSeq = seq.replace('T', 'U')


    if args.ntCounts and not args.gcContent and not args.rnaSeq and not args.fasta:
        print('Nucleotide count(A,C,G,T): ',ntCounts)
    if not args.ntCounts and args.gcContent and not args.rnaSeq and not args.fasta:
        print('gc content: ',gcContent)
    if not args.ntCounts and not args.gcContent and args.rnaSeq and not args.fasta:
        print('Transcription sequence: ',rnaSeq)

#若需要加print条目，则写在if里面(if/and not)
if not args.ntCounts and not args.gcContent and not args.rnaSeq and args.fasta:
    print('Sequence: ', fasta)

# if __name__ == '__main__':
    # for key in fasta.keys():
    #     print("name: ", key)
    #     print("sequence: ",fasta[key])
    #     print("nt_count: ", nt_count(fasta[key]))
    #     print("cg_content: ", cg_content(fasta[key]))
    #     print("rna: ", dna_trans_rna(fasta[key]))
    #     # print("protein: ", rna_trans_protein(fasta[key]))
    #     print("reverse_comple: ", reverse_comple("dna", fasta[key]))

def translate_dna(sequence):

    codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    proteinsequence = ''
    start = sequence.find('ATG')
    sequencestart = sequence[int(start):]
    stop = sequencestart.find('TAA')
    cds = str(sequencestart[:int(stop)+3])

#############################################################
#############################################################
#############################################################

import argparse
import os.path


def read_file(x):
    with open(x) as fa:
        sequences = {}
        for line in fa:
            if line.startswith(">"):
                name = line.rstrip("\n")
                sequences[name] = ""
            else:
                sequences[name] = sequences[name] + line.rstrip("\n")
    return sequences


def ntCounts_stat(fasta):
    for name in fasta.keys():
        ntCounts = []
        seq = fasta[name]
        for nt in ['A', 'C', 'G', 'T']:
            ntCounts.append(seq.count(nt))
    return ntCounts

def gcCount_stat(fasta):
    for name in fasta.keys():
        ntCounts = []
        seq = fasta[name]
    total = len(seq)
    gcCount = seq.count('G') + seq.count('C')
    gcContent = format(float(gcCount / total * 100), '.6f')
    return gcCount

def rna_seq(fasta):
    for name in fasta.keys():
        ntCounts = []
        seq = fasta[name]
    total = len(seq)
    gcCount = seq.count('G') + seq.count('C')
    rnaSeq = seq.replace('T', 'U')
    return rnaSeq

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-fa', action='store_true', dest='fasta', default=False, help='Check out the sequance')
    parser.add_argument('-nt', action='store_true', dest='ntCounts', default=False, help='Calculate Nucleotide count')
    parser.add_argument('-gc', action='store_true', dest='gcContent', default=False, help='Calculate gc content ')
    parser.add_argument('-r', action='store_true', dest='rnaSeq', default=False,
                        help='Calculate Transcription sequence')

    parser.add_argument('-i', dest='filename', type=str, help='Input file')  # 读取文件
    args = parser.parse_args()

    x = args.filename
    fasta = read_file(x)
    if args.ntCounts and not args.gcContent and not args.rnaSeq and not args.fasta:
        print('Nucleotide count(A,C,G,T): ',ntCounts_stat(fasta))
    if not args.ntCounts and args.gcContent and not args.rnaSeq and not args.fasta:
        print('gc content: ',gcCount_stat(fasta))
    if not args.ntCounts and not args.gcContent and args.rnaSeq and not args.fasta:
        print('Transcription sequence: ',rna_seq(fasta))

    if not args.ntCounts and not args.gcContent and not args.rnaSeq and args.fasta:
        print('Sequence: ', fasta)





