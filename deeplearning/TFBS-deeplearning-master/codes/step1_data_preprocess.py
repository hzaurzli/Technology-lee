import pandas as pd
import numpy as np
import tensorflow as tf
data=np.array(pd.read_csv('../data/knox6_train_data.csv'))
data_x=np.zeros((73520,200)).astype(np.str)
data_y=data[:,1]
for i in range(len(data[:,0])):
    data_x[i]=list(data[i,0])
data_x[data_x=='A']=3
data_x[data_x=='T']=2
data_x[data_x=='G']=1
data_x[data_x=='C']=0
data_x=data_x.astype(np.int)
data_x=np.transpose(tf.one_hot(data_x,depth=4),(0,2,1))
np.save('x.npy',data_x)
np.save('y.npy',data_y)
