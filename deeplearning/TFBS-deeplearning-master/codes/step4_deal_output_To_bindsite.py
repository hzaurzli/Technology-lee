import pandas as pd
import numpy as np
f1 = open("../result/knox6LSTM_bindingsite.txt",'w')
data1=np.array(pd.read_csv('../result/knox6LSTM.csv'))
data2=np.array(pd.read_csv('../data/greenleaf_ATAC_knox6predicted_seq.CSV'))
for i in range(len(data2[:,0])):
	chr1 = data2[i,0]
	start = data2[i,1]
	end = data2[i,2]
	if int(data1[i,0]) == 0:
		continue
	if int(data1[i,0]) != 0:
		start1 = int(start)+100*(int(data1[i,0])-1)
		end1 = int(start1) + 200
		f1.write(str(chr1)+"\t"+str(start1)+"\t"+str(end1)+"\n")
f1.close()
