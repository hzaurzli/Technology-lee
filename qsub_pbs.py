import os
import subprocess as sub

def qsub(pbsfile):
    cmd = "qsub " + pbsfile
    sub.call(cmd, shell=True)

# def qsub(ppn,pbsfile):
#     cmd = "qsub -q middle -V -l nodes=1:ppn="+str(ppn)+" "+pbsfile
#     sub.call(cmd, shell=True)

samples = []
for i in os.listdir("/home/rzli/swxxxjz/bam/10"):
    i = i.split("_")
    sample = i[0]
    samples.append(sample)
    samples = list(set(samples))

os.mkdir("/home/rzli/test/")
os.mkdir("/home/rzli/test/pbs")

for i in samples:
    with open("/home/rzli/test/pbs/" + i + ".sh", "w") as w:
        line = "#!/bin/bash" + "\n"
        line = "samtools view -s 0.5 -bo "
        line = line + "/home/rzli/swxxxjz/bam/10/" + i
        line = line + " /home/rzli/test/0.5_" + i
        w.write(line)

for i in samples:
    qsub("/home/rzli/test/pbs/" + i + ".sh")
