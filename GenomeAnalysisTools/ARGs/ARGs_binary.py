import sys,os


infile = open('/home/rzli/ssuis/LSSP/faa/data_results.out')
tmp = '/home/rzli/ssuis/LSSP/faa/tmp.txt'
out = '/home/rzli/ssuis/LSSP/faa/arg_binary.txt'

li_1 = []
li_2 = []
list_key = []
dic = {}
for i in infile:
    lis_1 = i.strip().split('\t')[1::]
    key = i.strip().split('\t')[0]
    for j in lis_1:
        j = j.split(':')
        li_1.append(j[0] + ':' + j[1])
        li_2.append(j[0] + ':' + j[1])
    dic[key] = li_2
    li_2 = []
    list_key.append(key)

list_1 = list(set(li_1))
print(dic['LSSP96'])

with open(tmp,'w') as w:
    for x in list_key:
        for y in list_1:
            line = x + '#' + y + '\n'
            w.write(line)
w.close()

f = open(tmp,'r')
with open(out,'w') as ww:
    for t in f:
        name = t.strip().split('#')[0]
        arg = t.strip().split('#')[1]
        if arg in dic[name]:
            line = name + '\t' + arg.strip().split('|')[2] + '\t' + '1' + '\n'
            ww.write(line)
            line = ''
        if arg not in dic[name]:
            line = name + '\t' + arg.strip().split('|')[2] + '\t' + '0' + '\n'
            ww.write(line)
            line = ''
ww.close()

os.remove(tmp)