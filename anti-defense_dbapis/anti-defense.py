#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  anti-defense.py 
#
#  Copyright 2025 Small runze
#  <small.runze@gmail.com> Small runze
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., HZAU.


import argparse
import os,sys,re,time
import random
import subprocess as sub
from subprocess import *
import subprocess as sub
import glob
import shutil
import biolib
from Bio import SeqIO
from Bio import AlignIO
from Bio import pairwise2 as pw2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis


class tools:
    def __init__(self):
        self.prokka = 'prokka'
        self.cdHit = 'cd-hit'
        self.hmmsearch = 'hmmsearch'


    def run(self, cmd, wkdir=None):
        sys.stderr.write("Running %s ...\n" % cmd)
        p = Popen(cmd, shell=True, cwd=wkdir)
        p.wait()
        return p.returncode

    def run_prokka(self, fastain, fastaout, prefix, type_annotation):
        cmd = '%s %s -o %s --prefix %s --kingdom %s -force' % (self.prokka, fastain, fastaout,prefix,type_annotation)
        return cmd

    def run_cdhit(self,inputfile, out, cutoff):
        cmd = '%s -i %s -o %s -c %s -M 0' % (self.cdHit, inputfile, out, cutoff)
        return cmd

    def run_hmmsearch(self,tblout, e_val, hmm, inputfile):
        cmd = '%s --tblout %s -E %s --cpu 2 %s %s' % (self.hmmsearch, tblout, e_val, hmm, inputfile)
        return cmd
        
    def run_hmmsearch_2(self,out, e_val, hmm, inputfile):
        cmd = '%s --domtblout %s -E %s --cpu 2 %s %s' % (self.hmmsearch, out, e_val, hmm, inputfile)
        return cmd

def molecular_weight(protein_fa,protein_filter_fa,MWU,MWL):
    protein_fa_info = open(protein_fa, "r")
    out_file = open(protein_filter_fa, "a")
    molecular_weight = open("./molecular_weight.txt", "w")
    for record in SeqIO.parse(protein_fa_info, "fasta"):
        ID_contig = record.id
        Seq_use = record.seq
        Desp = record.description
        protein_seq = str(Seq_use)
        if 'X' not in protein_seq and '*' not in protein_seq[:-1]:
            X = ProteinAnalysis(protein_seq)
            MW_cal = "%0.2f" % X.molecular_weight()
            if float(MW_cal) <= float(MWU) and float(MW_cal) >= float(MWL):
                element_record = SeqRecord(Seq_use, id=ID_contig, description=Desp)
                SeqIO.write(element_record, out_file, "fasta")
                molecular_weight.write(ID_contig + "\t" + MW_cal + "\n")

    protein_fa_info.close()
    molecular_weight.close()


def fasta2dict(fasta_name):
    with open(fasta_name) as fa:
        fa_dict = {}
        for line in fa:
            line = line.replace('\n', '')
            if line.startswith('>'):
              seq_name = line[0:]
              fa_dict[seq_name] = ''
            else:
              fa_dict[seq_name] += line.replace('\n', '')
    return fa_dict


def fasta2dict_2(fasta_name):
    with open(fasta_name) as fa:
        fa_dict = {}
        for line in fa:
            line = line.replace('\n', '')
            if line.startswith('>'):
              seq_name = line[1::].strip()
              fa_dict[seq_name] = ''
            else:
              fa_dict[seq_name] += line.replace('\n', '')
    return fa_dict


def find_pfam(cdhit_filter,lyase_list):
    input_file = open(cdhit_filter, "r")
    info_file = open('./hmmer_out/all_protein_filter_hmmer_out.txt', "r")
    info_file_two = open(lyase_list, "r")

    f2 = open(r'./hmmer_out/all_protein_filter_hmmer_out2.txt', 'w')
    for i in info_file:
        line = re.split('\s+', i)
        new_line = ' '.join(line)
        new_line = new_line.strip(' ')
        if new_line[0] != '#':
            f2.write(new_line)
            f2.write("\n")
    info_file.close()
    f2.close()

    info_pfam_list = []
    for i in info_file_two:
        data2 = i.strip()
        pfam_num_reported = str(data2)
        info_pfam_list.append(pfam_num_reported)

    info_file_new = open(r'./hmmer_out/all_protein_filter_hmmer_out2.txt', 'r')
    info_list = []
    for j in info_file_new:
        data1 = j.strip().split(' ')
        gen_id = str(data1[0])
        pfam_num = str(data1[2])
        for k in info_pfam_list:
            if k in pfam_num:
                info_list.append(gen_id)

    out_file = open("./all_protein_pfam_protein.fasta", "a")
    for record in SeqIO.parse(input_file, 'fasta'):
        Contig_ID = record.id
        Desp = record.description
        for i in info_list:
            gen_id = i
            if gen_id == Contig_ID:
                gene_seq = record.seq
                element_record = SeqRecord(gene_seq, id='', description=Desp)
                SeqIO.write(element_record, out_file, "fasta")
    input_file.close()
    info_file_new.close()
    info_file_two.close()
    out_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Lysin finder")
    parser.add_argument("-p", "--path", required=True, type=str, help="genome sequnce path")
    parser.add_argument("-t", "--type", required=True, type=str, help="prokka kingdom type")
    parser.add_argument("-pp", "--prophage_method", required=False, default='phispy', type=str, help="prophage predict method")
    parser.add_argument("-ds", "--dbscan_swa", required=False, type=str, default='', help="path of dbscan-swa.py")
    parser.add_argument("-c", "--cdhit_cutoff", default=0.95,required=False, type=float, help="cdhit cluster cutoff")
    parser.add_argument("-hc", "--hmmer_cutoff", default=1e-5,required=False, type=float, help="hmmer search cutoff")
    parser.add_argument("-hd", "--hmmer_db", required=True, type=str, help="reported lysin structures hmmer database path")
    parser.add_argument("-rl", "--reported_lysin", required=True, type=str, help="reported lysin structures(hmm files)")
    parser.add_argument("-wkdir", "--workdir", required=True, type=str, help="work directory")
    parser.add_argument("-mu", "--MWU", required=False, default=50000, type=float, help="upper proteins molecular weight")
    parser.add_argument("-ml", "--MWL", required=False, default=10000, type=float, help="lower proteins molecular weight")
    Args = parser.parse_args()
    
    
    if Args.workdir[-1] == '/':
        resultdir = os.path.basename(Args.workdir[:-1])
    elif Args.workdir[-1] == "\\":
        resultdir = os.path.basename(Args.workdir[:-1])
    else:
        resultdir = os.path.basename(Args.workdir)
    
    if os.path.isdir(os.path.dirname(os.path.abspath(Args.workdir)) +'/' + resultdir + '/') == True:
        pass
    else:
        os.mkdir(os.path.dirname(os.path.abspath(Args.workdir)) +'/' + resultdir + '/')
          

    tl = tools()        
    # step 1 prokka annotates ORFs
    curr_dir = sub.getoutput('pwd')
    os.chdir(Args.workdir)
    if os.path.isdir('./prokka_result/') == True:
        pass
    else:
        os.mkdir('./prokka_result/')

    target = Args.path
    curr_dir_target = curr_dir
    if target[-1] == '/':
        target = target
    elif target[-1] != '/':
        target = target + '/'

    if target[0] == '.':
        if target[1] == '/':
            target_suffix = target[1:]
        elif target[1] == '.':
            curr_dir_target = os.path.abspath(os.path.join(os.path.dirname(curr_dir + '/'), os.path.pardir))
            target_suffix = target[2:]
    else:
        target_suffix = target
        curr_dir_target = ''

    type_annotation = Args.type
    for i in os.listdir(curr_dir_target + target_suffix):
        lis = i.split('.')[:-1]
        name = '.'.join(lis)
        suffix = i.split('.')[-1]
        cmd_1 = tl.run_prokka(curr_dir_target + target_suffix + i,
                          './prokka_result/' + name + '/',name,type_annotation)
        tl.run(cmd_1)
        
    for i in os.listdir('./prokka_result/'):
      for j in os.listdir('./prokka_result/' + i):
        if j.endswith('.faa'):
          os.system('cat %s > %s' % ('./prokka_result/' + i + '/' + j, './prokka_result/' + i + '/tmp.txt'))
          with open('./prokka_result/' + i + '/' + j, 'w') as w:
            f = open('./prokka_result/' + i + '/tmp.txt')
            for line in f:
              if line.startswith('>'):
                print(line)
                first_line = line[1::].strip()
                key = first_line.split(' ')[0].split('_')[0]
                num = first_line.split(' ')[0].split('_')[1]
                lis = j.split('.')[:-1]
                name = '.'.join(lis)
                w.write('>' + name + '_' + num + '\n')
              else:
                w.write(line)
          w.close()

    

    # step 2 move faa into phage_faa fold
    if os.path.isdir('./phage_faa/') == True:
        pass
    else:
        os.mkdir('./phage_faa/')
        
    for i in os.listdir('./prokka_result/'):
        for j in os.listdir('./prokka_result/' + i):
            if os.path.splitext(j)[-1] == ".faa":
                os.system('cp %s %s' % ('./prokka_result/' + i + '/' + j, './phage_faa/'))

    
    if len(os.listdir('./phage_faa/')) == 0:
      os.system('rm -r ./prokka_result/ ./phage_faa/')
      raise('No phage faa found!')
      
    else:
      # step 3 phage faa together
      os.system('cat ./phage_faa/* > all_protein_ut.faa')
      
      fa_dict = fasta2dict('./all_protein_ut.faa')

      filters = ["B","Z","J","O","U","X",'*']
      with open('./all_protein.faa','w') as f:
          for key in fa_dict:
              if all(f not in fa_dict[key] for f in filters):
                  line = key + '\n' + fa_dict[key] + '\n'
                  f.write(line)
      f.close()
      
      if not os.path.getsize('./all_protein.faa'):
        with open('./putative_lysins.fa','w') as w:
          line = 'No lysins found!'
          w.write(line)
        w.close()
        
        os.system('rm -r ./prokka_result/')
        os.remove('./all_protein.faa')
        os.remove('./all_protein_ut.faa')
      
      else:
        # step 4 cdhit cluster
        cmd_4 = tl.run_cdhit('./all_protein.faa','./all_protein_cdhit.faa', Args.cdhit_cutoff)
        tl.run(cmd_4)

        # step 5 calculate molecular weight
        molecular_weight('./all_protein_cdhit.faa','./all_protein_cdhit_filter.faa', float(Args.MWU),float(Args.MWL))


        # step 6 hmmsearch reported lysin structure in pfam
        if os.path.isdir('./hmmer_out/') == True:
            pass
        else:
            os.mkdir('./hmmer_out/')

        hmmer_db = Args.hmmer_db
        if hmmer_db[0] == '.':
            if hmmer_db[1] == '/':
                hmmer_db_suffix = hmmer_db[1:]
                curr_dir_hmmerdb = curr_dir
            elif hmmer_db[1] == '.':
                curr_dir_hmmerdb = os.path.abspath(os.path.join(os.path.dirname(curr_dir + '/'), os.path.pardir))
                hmmer_db_suffix = hmmer_db[2:]
        else:
            hmmer_db_suffix = hmmer_db
            curr_dir_hmmerdb = ''

        cmd_5 = tl.run_hmmsearch('./hmmer_out/all_protein_filter_hmmer_out.txt', Args.hmmer_cutoff,
                                 curr_dir_hmmerdb + hmmer_db_suffix,
                                 './all_protein_cdhit_filter.faa')
        tl.run(cmd_5)
        
        cmd_5_p = tl.run_hmmsearch_2('./hmmer_out/all_protein.txt', Args.hmmer_cutoff,
                                   curr_dir_hmmerdb + hmmer_db_suffix,
                                   './all_protein_cdhit_filter.faa')
        tl.run(cmd_5_p)

        reported_lysin = Args.reported_lysin
        if reported_lysin[0] == '.':
            if reported_lysin[1] == '/':
                reported_lysin_suffix = reported_lysin[1:]
                curr_dir_rp = curr_dir
            elif reported_lysin[1] == '.':
                curr_dir_rp = os.path.abspath(os.path.join(os.path.dirname(curr_dir + '/'), os.path.pardir))
                reported_lysin_suffix = reported_lysin[2:]
        else:
            reported_lysin_suffix = reported_lysin
            curr_dir_rp = ''

        find_pfam('./all_protein_cdhit_filter.faa', curr_dir_rp + reported_lysin_suffix)
        

        # step 9 combine results of CAZY and pfam
        os.system('cat all_protein_pfam_protein.fasta > pfam_anti-defense.fasta')
        
        dic_fa = fasta2dict_2('./pfam_anti-defense.fasta')

        f1 = open('./molecular_weight.txt')
        f2 = open('./hmmer_out/all_protein.txt')
        with open('./MW_Length.txt', 'w') as w1:
          for i in f1:
            name = i.strip().split('\t')[0]
            mw = i.strip().split('\t')[1]
            if name in dic_fa.keys():
              line = name + '\t' + mw + '\t' + str(len(dic_fa[name])) + '\n'
              w1.write(line)
        w1.close()
        

        Domain_Info_lis = []
        with open('./Domain_Info.txt', 'w') as w2:
          for line in f2:
            if line[0] != "#" and len(line.split())!=0:
              arr = line.strip().split(" ")
              arr = list(filter(None, arr))
              name = arr[0]
              if name in dic_fa.keys():
                li = arr[0] + '\t' + arr[3] + '(Length:' + arr[5] + ')' + '\t' + arr[21] + '\t' + arr[19] + '-' + arr[20] + '\n'
                print(li)
                Domain_Info_lis.append(li)
         
          Domain_Info_lis_new = list(set(Domain_Info_lis))
          for line in Domain_Info_lis_new:
            w2.write(line)
        w2.close()
        
        # os.system("sed -i '$d' %s" % ('/home/runzeli/rzli/zy/result/Domain_Info.txt'))
        
        f1 = open('./MW_Length.txt')
        f2 = open('./Domain_Info.txt')
        
        
        dic_info = {}
        for lines in f1:
          line = lines.strip().split('\t')
          id_1 = line[0]
          mw = line[1]
          length = line[2]
          mw_length = []
          mw_length.append(mw)
          mw_length.append(length)
          dic_info[id_1] = mw_length
          
        
        for lines in f2:
          line = lines.strip().split('\t')
          id_2 = line[0]
          pf = line[1] + '&' + line[2] + '&' + line[3]
          if id_2 in dic_info.keys():
            dic_info[id_2].append(pf)
        
        print(dic_info)
        os.system('cat pfam_anti-defense.fasta > anti-defense.fasta')
        with open('./anti-defense_info.txt','w') as w:
            line = 'ID' + '\t' + 'MW' + '\t' + 'Length' + '\t' + 'Domains' + '\n'
            w.write(line)
            for key in dic_info:
              line = key + '\t' + '\t'.join(dic_info[key][0:2]) + '\t' + ';'.join(dic_info[key][2:len(dic_info[key])]) + '\n'
              w.write(line)
        w.close()
              
              

