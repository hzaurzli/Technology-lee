import sys,os


infile = '/home/rzli/ssuis/card/index/protein_ARGs.fasta'
outfile = '/home/rzli/ssuis/card/index/protein_ARGs_length.fasta'
sequence = ' '
fasta = {}

with open(infile) as file_one:
    for line in file_one:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            active_sequence_name = line[1:].split(' ')[0]
            if active_sequence_name not in fasta:
                fasta[active_sequence_name] = []
            continue
        sequence = line
        fasta[active_sequence_name].append(sequence)

with open(outfile,'w') as w:
    for key in fasta:
        print(len(fasta[key][0]))
        line = '>' + key + ':' + str(len(fasta[key][0])) + '\n' + fasta[key][0] + '\n'
        w.write(line)
