#同步
from multiprocessing import Pool
import time
def test(p):
       print(p)
       time.sleep(3)
if __name__=="__main__":
    pool = Pool(processes=10)
    for i in range(500):
        '''
        ('\n'
         '	（1）遍历500个可迭代对象，往进程池放一个子进程\n'
         '	（2）执行这个子进程，等子进程执行完毕，再往进程池放一个子进程，再执行。（同时只执行一个子进程）\n'
         '	 for循环执行完毕，再执行print函数。\n'
         '	')
        '''
        pool.apply(test, args=(i,))
    print('test')
    pool.close()
    pool.join()

#异步（非阻塞）
from multiprocessing import Pool
import time
def test(p):
       print(p)
       time.sleep(3)
if __name__=="__main__":
    pool = Pool(processes=2)
    for i in range(500):
        '''
         （1）循环遍历，将500个子进程添加到进程池（相对父进程会阻塞）\n'
         （2）每次执行2个子进程，等一个子进程执行完后，立马启动新的子进程。（相对父进程不阻塞）\n'
        '''
        pool.apply_async(test, args=(i,))   #维持执行的进程总数为10，当一个进程执行完后启动一个新进程.
    print('test')
    pool.close()
    pool.join()

#map
from multiprocessing import Pool
def test(i):
    print(i)
if  __name__ == "__main__":
    lists = [1, 2, 3]
    pool = Pool(processes=2)       #定义最大的进程数
    pool.map(test, lists)          #p必须是一个可迭代变量。
    pool.close()
    pool.join()


##########
from multiprocessing import Pool
import time
import os , sys

samples = []
for i in os.listdir("/home/rzli/swxxxjz/bam/10"):
    i = i.split("_")
    sample = i[0]
    samples.append(sample)
    samples = list(set(samples))

def test(a,b):
    os.system("samtools view -s 0.5 -bo %s %s" % (a,b)) #占位符


if __name__=="__main__":
    pool = Pool(processes=5)
    for p in samples:
        a = "/home/rzli/test/0.5_" + p
        b = "/home/rzli/swxxxjz/bam/10/" + p
        pool.apply_async(test, args=(a,b))
      print('test')
    pool.close()
    pool.join()
