import os
import sys
import re
import argparse


def remove_fa(infile,outfile):
        with open(infile) as f:
            Dict = {}
            for line in f:
                if line[0] == ">":
                    key = line.strip()
                    Dict[key] = []
                else:
                    Dict[key].append(line.strip())

        with open(outfile, 'w') as o:
            for key, value in Dict.items():
                o.write("{}\n{}\n".format(key, ''.join(value)))


def fa_dict(input):
        f=open(input)
        seq={}
        for line in f:
                if line.startswith('>'):
                        name=line.replace('>','').split()[0]
                        seq[name]=''
                else:
                        seq[name]+=line.replace('\n','').strip()
        f.close()
        return seq


def main(infile,out):
        outfile = os.path.basename(infile).split('.')[0]
        remove_fa(infile, outfile+'tmp.fasta')
        dict = fa_dict(outfile+'tmp.fasta')
        os.remove(outfile+'tmp.fasta')

        lis = []
        for key in dict:
                if len(dict[key]) > 300:
                        lis.append(dict[key])

        with open(out,'w') as w:
                path = os.path.basename(out)
                name ='>' + path.split('.')[0] + '\n'
                seq = ''.join(lis) + '\n'
                line = name + seq
                w.write(line)
        w.close()


if __name__ == "__main__":
        parser = argparse.ArgumentParser(description="phage assembly")
        parser.add_argument("-i", "--input", required=True, type=str, help="")
        parser.add_argument("-o", "--output", required=True, type=str, help="")
        Args = parser.parse_args()

        main(Args.input,Args.output)