import argparse
import os,sys,re
import subprocess as sub
from subprocess import *


class tools:
    def __init__(self):
        self.makedb = 'makeblastdb'
        self.blastn = 'blastn'
        self.blastp = 'blastp'

        self.bedtools = 'bedtools'
        self.seqkit = 'seqkit'

    def run(self, cmd, wkdir=None):
        sys.stderr.write("Running %s ...\n" % cmd)
        p = Popen(cmd, shell=True, cwd=wkdir)
        p.wait()
        return p.returncode

    def run_blastn_makedb(self,fastain,fastaout):
        cmd = '%s -in %s -dbtype nucl -out %s' % (self.makedb, fastain, fastaout)
        return cmd


    def run_blastp_makedb(self,fastain,fastaout):
        cmd = '%s -in %s -dbtype prot -out %s' % (self.makedb, fastain, fastaout)
        return cmd

    def run_blastn(self,db,fastain,fastaout):
        cmd = '%s -db %s -query %s -out %s -outfmt 6 -evalue 1e-5' % (self.blastn, db, fastain, fastaout)
        return cmd

    def run_blastp(self, db,fastain,fastaout):
        cmd = '%s -db %s -query %s -out %s -outfmt 6 -evalue 1e-5' % (self.blastp, db, fastain, fastaout)
        return cmd


    def run_bedtools(self, a_bed, b_bed):
        cmd = '%s intersect -a %s -b %s -wo > overlap.txt' % (self.bedtools, a_bed, b_bed)
        return cmd

    def run_seqkit(self, fasta):
        cmd = '%s fx2tab --length --name --header-line %s -o ./dbFasta_length.txt' % (self.seqkit,fasta)
        return cmd

def main():
    parser = argparse.ArgumentParser(description="Blast")
    parser.add_argument("-i",
                        "--input",
                        required=True,
                        type=str,
                        help="All target genome")
    parser.add_argument("-d",
                        "--db_fa",
                        required=True,
                        type=str,
                        help="Blast db fasta")
    parser.add_argument("-t",
                        "--type",
                        required=True,
                        type=str,
                        help="Sequence type")
    parser.add_argument("-p",
                        "--prefix",
                        required=True,
                        type=str,
                        help="Blast db prefix")
    parser.add_argument("-s",
                        "--score",
                        required=True,
                        type=float,
                        help="Similarity score")
    parser.add_argument("-c",
                        "--coverage",
                        required=True,
                        type=float,
                        help="Coverage score")
    parser.add_argument("-l",
                        "--overlap_cov",
                        required=True,
                        type=float,
                        help="Overlap coverage percentage")
    parser.add_argument("-wkdir",
                        "--workdir",
                        required=True,
                        type=str,
                        help="Work directory")
    Args = parser.parse_args()


    tl = tools()

    curr_dir = sub.getoutput('pwd')
    curr_dir_input = sub.getoutput('pwd')
    curr_dir_db = sub.getoutput('pwd')

    target = Args.workdir
    if target[-1] == '/':
        target = target
    elif target[-1] != '/':
        target = target + '/'

    if target[0] == '.':
        if target[1] == '/':
            target_suffix = target[1:]
        elif target[1] == '.':
            curr_dir = os.path.abspath(os.path.join(os.path.dirname(curr_dir + '/'), os.path.pardir))
            target_suffix = target[2:]
    else:
        target_suffix = target
        curr_dir = ''

    if os.path.isdir(curr_dir + target_suffix) == True:
        pass
    else:
        os.mkdir(curr_dir + target_suffix)

    os.chdir(curr_dir + target_suffix)

    ##### count fa length

    target_input = Args.input
    if target_input[-1] == '/':
        target_input = target_input
    elif target_input[-1] != '/':
        target_input = target_input + '/'

    if target_input[0] == '.':
        if target_input[1] == '/':
            target_input_suffix = target_input[1:]
        elif target_input[1] == '.':
            curr_dir_input = os.path.abspath(os.path.join(os.path.dirname(curr_dir_input + '/'), os.path.pardir))
            target_input_suffix = target_input[2:]
    else:
        target_input_suffix = target_input
        curr_dir_input = ''

    dbFasta = open(Args.db_fa)
    cmd = tl.run_seqkit(Args.db_fa)
    tl.run(cmd)

    fr = open('./dbFasta_length.txt')
    dic = {}
    keys = []
    for line in fr:
        if line.startswith('#'):
            pass
        else:
            v = line.strip().split('\t')
            m = v[0].split(' ')
            dic[m[0]] = v[1]
            keys.append(m[0])
    fr.close()

    ### add the length into fasta
    with open('./length_dbFasta.fa','w') as w:
        for i in dbFasta:
            if i.startswith('>'):
                j = i.strip().split()[0][1:]
                j = '>' + j + '::' + dic[j] + '\n'
                w.write(j)
            else:
                j = i
                w.write(j)
    w.close()

    ##########length_fa
    db_fa = './length_dbFasta.fa'
    ### make blastdb
    if Args.type == 'n':
        if os.path.isdir('./' + Args.prefix) == True:
            pass
        else:
            os.mkdir('./' + Args.prefix)

        cmd = tl.run_blastn_makedb(db_fa,'./' + Args.prefix + '/' + Args.prefix)
        tl.run(cmd)

        cmd = tl.run_blastn('./' + Args.prefix + '/' + Args.prefix,
                            curr_dir_input + target_input_suffix[:-1],
                            './' + Args.prefix + '.m6')
        tl.run(cmd)

        blast_result = open('./' + Args.prefix + '.m6')
        with open('./' + Args.prefix + '_filter.m6','w') as w:
            for i in blast_result:
                length = int(i.strip().split('\t')[1].split('::')[1])
                length_cov = abs(int(i.strip().split('\t')[7]) - int(i.strip().split('\t')[6])) / length
                similarity = i.strip().split('\t')[2]
                if float(length_cov) > float(Args.coverage) and float(similarity)/100 > float(Args.score):
                    lines = i.strip().split('\t')[0] + '\t' + i.strip().split('\t')[6] + '\t' + \
                            i.strip().split('\t')[7] + '\t' + i.strip().split('\t')[2] + '\t' + \
                            i.strip().split('\t')[1] + '\n'
                    w.write(lines)
        w.close()

        cmd = tl.run_bedtools('./' + Args.prefix + '_filter.m6','./' + Args.prefix + '_filter.m6')
        tl.run(cmd)

        overlap_result = open('./overlap.txt')
        with open('./overlap_f.txt', 'w') as w:
            for i in overlap_result:
                name_1 = i.strip().split('\t')[0] + '-' + i.strip().split('\t')[1] + '-' + i.strip().split('\t')[2]
                name_2 = i.strip().split('\t')[5] + '-' + i.strip().split('\t')[6] + '-' + i.strip().split('\t')[7]
                item_1 = i.strip().split('\t')[4]
                item_2 = i.strip().split('\t')[9]

                if name_1 == name_2:
                    item = i.strip().split('\t')[4]
                    name = name_1
                    score = float(i.strip().split('\t')[3])

                    lines = name + '\t' + item + '\t' + str(score) + '\n'

                elif name_1 != name_2 and item_1 == item_2:
                    overlap_len = int(i.strip().split('\t')[10])
                    gene_len_1 = abs(int(i.strip().split('\t')[2]) - int(i.strip().split('\t')[1]))
                    gene_len_2 = abs(int(i.strip().split('\t')[7]) - int(i.strip().split('\t')[6]))
                    over_1 = float(overlap_len / gene_len_1)
                    over_2 = float(overlap_len / gene_len_2)
                    score_1 = float(i.strip().split('\t')[3])
                    score_2 = float(i.strip().split('\t')[8])

                    if over_1 < float(Args.overlap_cov) and over_2 < float(Args.overlap_cov):
                        lines = name_1 + '\t' + item_1 + '\t' + str(score_1) + '\n' \
                                + name_2 + '\t' + item_2 + '\t' + str(score_2) + '\n'

                    elif over_1 > float(Args.overlap_cov) and over_2 < float(Args.overlap_cov):
                        if score_1 > score_2 or score_1 == score_2:
                            item = i.strip().split('\t')[4]
                            name = name_1
                            score = score_1
                        else:
                            item = i.strip().split('\t')[9]
                            name = name_2
                            score = score_2

                        lines = name + '\t' + item + '\t' + str(score) + '\n'

                    elif over_1 < float(Args.overlap_cov) and over_2 > float(Args.overlap_cov):
                        if score_1 > score_2 or score_1 == score_2:
                            item = i.strip().split('\t')[4]
                            name = name_1
                            score = score_1
                        else:
                            item = i.strip().split('\t')[9]
                            name = name_2
                            score = score_2

                        lines = name + '\t' + item + '\t' + str(score) + '\n'

                    elif over_1 > float(Args.overlap_cov) and over_2 > float(Args.overlap_cov):
                        if score_1 > score_2 or score_1 == score_2:
                            item = i.strip().split('\t')[4]
                            name = name_1
                            score = score_1
                        else:
                            item = i.strip().split('\t')[9]
                            name = name_2
                            score = score_2

                        lines = name + '\t' + item + '\t' + str(score) + '\n'


                elif name_1 != name_2 and item_1 != item_2:
                    overlap_len = int(i.strip().split('\t')[10])
                    gene_len_1 = abs(int(i.strip().split('\t')[2]) - int(i.strip().split('\t')[1]))
                    gene_len_2 = abs(int(i.strip().split('\t')[7]) - int(i.strip().split('\t')[6]))
                    over_1 = float(overlap_len / gene_len_1)
                    over_2 = float(overlap_len / gene_len_2)
                    score_1 = float(i.strip().split('\t')[3])
                    score_2 = float(i.strip().split('\t')[8])

                    if over_1 < float(Args.overlap_cov) and over_2 < float(Args.overlap_cov):
                        lines = name_1 + '\t' + item_1 + '\t' + str(
                            score_1) + '\n' + name_2 + '\t' + item_2 + '\t' + str(score_2) + '\n'
                    elif over_1 > float(Args.overlap_cov) and over_2 < float(Args.overlap_cov):
                        if score_1 > score_2 or score_1 == score_2:
                            item = i.strip().split('\t')[4]
                            name = name_1
                            score = score_1
                        else:
                            item = i.strip().split('\t')[9]
                            name = name_2
                            score = score_2

                        lines = name + '\t' + item + '\t' + str(score) + '\n'

                    elif over_1 < float(Args.overlap_cov) and over_2 > float(Args.overlap_cov):
                        if score_1 > score_2 or score_1 == score_2:
                            item = i.strip().split('\t')[4]
                            name = name_1
                            score = score_1
                        else:
                            item = i.strip().split('\t')[9]
                            name = name_2
                            score = score_2

                        lines = name + '\t' + item + '\t' + str(score) + '\n'

                    elif over_1 > float(Args.overlap_cov) and over_2 > float(Args.overlap_cov):
                        if score_1 > score_2 or score_1 == score_2:
                            item = i.strip().split('\t')[4]
                            name = name_1
                            score = score_1
                        else:
                            item = i.strip().split('\t')[9]
                            name = name_2
                            score = score_2

                        lines = name + '\t' + item + '\t' + str(score) + '\n'

                w.write(lines)
        w.close()

        os.system('cat %s | sort | uniq > %s' % ('./overlap_f.txt', './overlap_filter.txt'))
        os.system('rm %s %s %s %s' % ('./overlap_f.txt', './overlap.txt',
                                      './dbFasta_length.txt', './length_dbFasta.fa'))


    elif Args.type == 'p':
        if os.path.isdir('./' + Args.prefix) == True:
            pass
        else:
            os.mkdir('./' + Args.prefix)

        cmd = tl.run_blastp_makedb(db_fa,'./' + Args.prefix + '/' + Args.prefix)
        tl.run(cmd)

        cmd = tl.run_blastp('./' + Args.prefix + '/' + Args.prefix,
                            curr_dir_input + target_input_suffix[:-1],
                            './' + Args.prefix + '.m6')
        tl.run(cmd)

        blast_result = open('./' + Args.prefix + '.m6')
        with open('./' + Args.prefix + '_filter.m6','w') as w:
            for i in blast_result:
                length = int(i.strip().split('\t')[1].split('::')[1])
                length_cov = abs(int(i.strip().split('\t')[7]) - int(i.strip().split('\t')[6])) / length
                similarity = i.strip().split('\t')[2]
                if float(length_cov) > float(Args.coverage) and float(similarity)/100 > float(Args.score):
                    lines = i.strip().split('\t')[0] + '\t' + i.strip().split('\t')[6] + '\t' + \
                            i.strip().split('\t')[7] + '\t' + i.strip().split('\t')[2] + '\t' + \
                            i.strip().split('\t')[1] + '\n'
                    w.write(lines)
        w.close()

        cmd = tl.run_bedtools('./' + Args.prefix + '_filter.m6','./' + Args.prefix + '_filter.m6')
        tl.run(cmd)

        overlap_result = open('./overlap.txt')
        with open('./overlap_f.txt', 'w') as w:
            for i in overlap_result:
                name_1 = i.strip().split('\t')[0] + '-' + i.strip().split('\t')[1] + '-' + i.strip().split('\t')[2]
                name_2 = i.strip().split('\t')[5] + '-' + i.strip().split('\t')[6] + '-' + i.strip().split('\t')[7]
                item_1 = i.strip().split('\t')[4]
                item_2 = i.strip().split('\t')[9]

                if name_1 == name_2:
                    item = i.strip().split('\t')[4]
                    name = name_1
                    score = float(i.strip().split('\t')[3])

                    lines = name + '\t' + item + '\t' + str(score) + '\n'

                elif name_1 != name_2 and item_1 == item_2:
                    overlap_len = int(i.strip().split('\t')[10])
                    gene_len_1 = abs(int(i.strip().split('\t')[2]) - int(i.strip().split('\t')[1]))
                    gene_len_2 = abs(int(i.strip().split('\t')[7]) - int(i.strip().split('\t')[6]))
                    over_1 = float(overlap_len / gene_len_1)
                    over_2 = float(overlap_len / gene_len_2)
                    score_1 = float(i.strip().split('\t')[3])
                    score_2 = float(i.strip().split('\t')[8])

                    if over_1 < float(Args.overlap_cov) and over_2 <  float(Args.overlap_cov):
                        lines = name_1 + '\t' + item_1 + '\t' + str(score_1) + '\n' \
                                + name_2 + '\t' + item_2 + '\t' + str(score_2) + '\n'

                    elif over_1 > float(Args.overlap_cov) and over_2 < float(Args.overlap_cov):
                        if score_1 > score_2 or score_1 == score_2:
                            item = i.strip().split('\t')[4]
                            name = name_1
                            score = score_1
                        else:
                            item = i.strip().split('\t')[9]
                            name = name_2
                            score = score_2

                        lines = name + '\t' + item + '\t' + str(score) + '\n'

                    elif over_1 < float(Args.overlap_cov) and over_2 > float(Args.overlap_cov):
                        if score_1 > score_2 or score_1 == score_2:
                            item = i.strip().split('\t')[4]
                            name = name_1
                            score = score_1
                        else:
                            item = i.strip().split('\t')[9]
                            name = name_2
                            score = score_2

                        lines = name + '\t' + item + '\t' + str(score) + '\n'

                    elif over_1 > float(Args.overlap_cov) and over_2 > float(Args.overlap_cov):
                        if score_1 > score_2 or score_1 == score_2:
                            item = i.strip().split('\t')[4]
                            name = name_1
                            score = score_1
                        else:
                            item = i.strip().split('\t')[9]
                            name = name_2
                            score = score_2

                        lines = name + '\t' + item + '\t' + str(score) + '\n'


                elif name_1 != name_2 and item_1 != item_2:
                    overlap_len = int(i.strip().split('\t')[10])
                    gene_len_1 = abs(int(i.strip().split('\t')[2]) - int(i.strip().split('\t')[1]))
                    gene_len_2 = abs(int(i.strip().split('\t')[7]) - int(i.strip().split('\t')[6]))
                    over_1 = float(overlap_len / gene_len_1)
                    over_2 = float(overlap_len / gene_len_2)
                    score_1 = float(i.strip().split('\t')[3])
                    score_2 = float(i.strip().split('\t')[8])

                    if over_1 < float(Args.overlap_cov) and over_2 < float(Args.overlap_cov):
                        lines = name_1 + '\t' + item_1 + '\t' + str(
                            score_1) + '\n' + name_2 + '\t' + item_2 + '\t' + str(score_2) + '\n'
                    elif over_1 > float(Args.overlap_cov) and over_2 < float(Args.overlap_cov):
                        if score_1 > score_2 or score_1 == score_2:
                            item = i.strip().split('\t')[4]
                            name = name_1
                            score = score_1
                        else:
                            item = i.strip().split('\t')[9]
                            name = name_2
                            score = score_2

                        lines = name + '\t' + item + '\t' + str(score) + '\n'

                    elif over_1 < float(Args.overlap_cov) and over_2 > float(Args.overlap_cov):
                        if score_1 > score_2 or score_1 == score_2:
                            item = i.strip().split('\t')[4]
                            name = name_1
                            score = score_1
                        else:
                            item = i.strip().split('\t')[9]
                            name = name_2
                            score = score_2

                        lines = name + '\t' + item + '\t' + str(score) + '\n'

                    elif over_1 > float(Args.overlap_cov) and over_2 > float(Args.overlap_cov):
                        if score_1 > score_2 or score_1 == score_2:
                            item = i.strip().split('\t')[4]
                            name = name_1
                            score = score_1
                        else:
                            item = i.strip().split('\t')[9]
                            name = name_2
                            score = score_2

                        lines = name + '\t' + item + '\t' + str(score) + '\n'

                w.write(lines)
        w.close()

        os.system('cat %s | sort | uniq > %s' % ('./overlap_f.txt', './overlap_filter.txt'))
        os.system('rm %s %s %s %s' % ('./overlap_f.txt', './overlap.txt',
                                './dbFasta_length.txt','./length_dbFasta.fa'))


if __name__ == "__main__":
    main()



