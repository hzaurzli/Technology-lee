import joblib
import pandas as pd
import numpy as np
import tensorflow as tf
from keras.models import load_model
data=np.array(pd.read_csv('../data/greenleaf_ATAC_knox6predicted_seq.csv'))
def pred(model,data):
    result=np.zeros(data.shape[0])
    for i in range(data.shape[0]):
        seq = list(data[i, 3])
        length = len(seq)
        seqarr = np.zeros((1, 200)).astype(np.str)
        for a in range(length // 200 - 1):
            seqarr[0] = seq[a * 200:(a + 1) * 200]
            seqarr[ seqarr== 'A'] = 3
            seqarr[ seqarr == 'T'] = 2
            seqarr[ seqarr == 'G'] = 1
            seqarr[ seqarr == 'C'] = 0
            seq2arr = np.array(seqarr).astype(np.int)
            datax=np.transpose(tf.one_hot(seq2arr,depth=4),(0,2,1))
            type=model.predict(datax)
            if type>=0.5:
                result[i]=a+1
                break
            else:
                result[i]=0
    return result
def svmpred(model,data):
    result=np.zeros(data.shape[0])
    for i in range(data.shape[0]):
        seq = list(data[i, 3])
        length = len(seq)
        seqarr = np.zeros((1, 200)).astype(np.str)
        for a in range(length // 200 - 1):
            seqarr[0] = seq[a * 200:(a + 1) * 200]
            seqarr[ seqarr== 'A'] = 3
            seqarr[ seqarr == 'T'] = 2
            seqarr[ seqarr == 'G'] = 1
            seqarr[ seqarr == 'C'] = 0
            seq2arr = np.array(seqarr).astype(np.int)
            datax=np.transpose(tf.one_hot(seq2arr,depth=4),(0,2,1))
            datax=datax.reshape(-1,4*200)
            type=model.predict(datax)
            if type>=0.5:
                result[i]=a+1
                break
            else:
                result[i]=0
    return result

TCNmodel = load_model('knox6TCN.h5')
TCNresult = pd.DataFrame(pred(TCNmodel,data)).to_csv('knox6TCN.csv')
LSTMmodel = load_model('knox6LSTM.h5')
LSTMresult = pd.DataFrame(pred(LSTMmodel,data)).to_csv('knox6LSTM.csv')
SVMmodel=joblib.load('knox6svm.h5')
svmresult=pd.DataFrame(svmpred(SVMmodel,data)).to_csv('knox6SVM.csv')