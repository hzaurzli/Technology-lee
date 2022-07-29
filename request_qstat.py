import os
import subprocess as sub
import re
import time

samples = []
for i in os.listdir("/home/rzli/swxxxjz/codes"):
    samples.append(i)
    samples = list(set(samples))

def qsub(pbsfile):
    cmd = "qsub " + pbsfile
    sub.call(cmd, shell=True)

def qstat(x = 1):
    c = sub.getoutput('qstat')
    a = re.compile(r'admin1')
    b = a.findall(c)
    if len(b) != 0:
        time.sleep(10)
        return qstat(x=1)
    else:
        return 1

if __name__ == "__main__":
    for i in samples[0:6]:
        qsub("/home/rzli/swxxxjz/codes/" + i)
    if qstat(x = 1) == 1:
        for i in samples[6:12]:
            qsub("/home/rzli/swxxxjz/codes/" + i)
