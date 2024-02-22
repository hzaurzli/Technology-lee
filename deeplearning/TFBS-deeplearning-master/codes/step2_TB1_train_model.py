import os
import  matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
from plotnine import *
import keras.optimizers
from tensorflow.keras.models import Model
from tensorflow.keras.layers import add, Input, Conv1D, Embedding, Flatten, Dense,BatchNormalization,MaxPooling1D,Dropout,Bidirectional,LSTM,GlobalAveragePooling1D
from tensorflow.keras import regularizers
from keras.layers import ReLU
from tensorflow.keras import initializers
from sklearn.svm import SVC,LinearSVC
from tensorflow.keras.layers import Multiply
from tensorflow.python.keras.layers.core import *
from sklearn.metrics import accuracy_score
import tensorflow.keras.backend as K
from sklearn.model_selection import train_test_split
from tensorflow.keras.models import *
from tensorflow.keras.metrics import Accuracy
from tensorflow.keras.losses import binary_crossentropy
import joblib
# from thundersvm import SVC
def train_svm():
    x = np.load('x.npy').astype(np.float16)
    x=x.reshape(-1,200*4)
    y = np.load('y.npy', allow_pickle=True).astype(np.int)
    train_x, test_x, train_y, test_y = train_test_split(x, y, test_size=0.2, random_state=43)
    model=SVC(probability=True,kernel='linear',verbose=1)
    model.fit(train_x,train_y)
    score=accuracy_score(model.predict(test_x),test_y)
    joblib.dump(model,'svm.h5')

def ResBlock(x, filters, kernel_size, dilation_rate):
    r = Conv1D(filters, kernel_size, padding='same', dilation_rate=dilation_rate,kernel_initializer=initializers.RandomNormal(),kernel_regularizer=regularizers.l2(0.0002))(x)  # 第一卷积
    r = BatchNormalization()(r)
    r = ReLU()(r)
    r = Conv1D(filters, kernel_size, padding='same', dilation_rate=dilation_rate,kernel_initializer=initializers.RandomNormal())(r)  # 第二卷积
    r = BatchNormalization()(r)
    r = ReLU()(r)
    if x.shape[-1] == filters:
        shortcut = x
    else:
        shortcut = Conv1D(filters, kernel_size, padding='same')(x)  # shortcut（捷径）
    o = add([r, shortcut])
    o = Activation('relu')(o)  # **函数
    return o
def TCN():
    inputs = Input(shape=(4,200))
    x = ResBlock(inputs, filters=16, kernel_size=3, dilation_rate=1)###TCN
    x = ResBlock(x, filters=32, kernel_size=3, dilation_rate=3)###TCN
    x = ResBlock(x, filters=64, kernel_size=3, dilation_rate=9)###TCN
    x = ResBlock(x, filters=128, kernel_size=3, dilation_rate=27)###TCN
    x = GlobalAveragePooling1D()(x)
    x = Dense(50,activation='relu',kernel_initializer=initializers.RandomNormal())(x)
    x = Dropout(0.1)(x)
    x = Dense(1,activation='sigmoid')(x)
    model = Model(inputs, x)
    return model
def LSTMs():
    inputs = Input(shape=(4, 200))
    x = Bidirectional(LSTM(50,return_sequences=True,activation='tanh'))(inputs) ###TCN
    x = Dropout(0.02)(x)
    x = Bidirectional(LSTM(50,return_sequences=False,activation='tanh'))(x)  ###TCN
    x = Dropout(0.02)(x)
    x = Bidirectional(LSTM(50,return_sequences=False,activation='tanh'))(x)  ###TCN
    x = Dropout(0.02)(x)
    x = Dense(50, activation='relu', kernel_initializer=initializers.RandomNormal())(x)
    x = Dropout(0.1)(x)
    x = Dense(1, activation='sigmoid')(x)
    model = Model(inputs, x)
    return model
class TCN_Net():
    def train(self, epochs, batch_size=128, sample_interval=50):
        x=np.load('x.npy').astype(np.float)
        y=np.load('y.npy',allow_pickle=True).astype(np.int)
        train_x,test_x,train_y,test_y=train_test_split(x,y,test_size=0.2,random_state=40)
        TCNmodel=TCN()
        TCNmodel.compile(optimizer='adam', loss=binary_crossentropy, metrics=['accuracy'])
        TCNmodel.summary()
        from keras.callbacks import ModelCheckpoint
        filepath = 'TCN.h5'  ####原最佳模型7dimmodelATT_TCNnew.h5
        checkpoint = ModelCheckpoint(filepath,
                                     monitor='val_loss',
                                     verbose=0,
                                     save_best_only=True,
                                     save_weights_only=False,
                                     mode='auto',
                                     period=1
                                     )
        callbacks_list = [checkpoint]
        cost = TCNmodel.fit(train_x, train_y, batch_size=batch_size, epochs=epochs
                         , verbose=1, validation_data=(test_x, test_y), callbacks=callbacks_list)
class LSTM_Net():
    def train(self, epochs, batch_size=128, sample_interval=50):
        x=np.load('x.npy').astype(np.float)
        y=np.load('y.npy',allow_pickle=True).astype(np.int)
        train_x,test_x,train_y,test_y=train_test_split(x,y,test_size=0.2,random_state=10)
        LSTMmodel=LSTMs()
        LSTMmodel.compile(optimizer='adam', loss=binary_crossentropy, metrics=['accuracy'])
        LSTMmodel.summary()
        from keras.callbacks import ModelCheckpoint
        filepath = 'LSTM.h5'  ####原最佳模型7dimmodelATT_TCNnew.h5
        checkpoint = ModelCheckpoint(filepath,
                                     monitor='val_loss',
                                     verbose=0,
                                     save_best_only=True,
                                     save_weights_only=False,
                                     mode='auto',
                                     period=1
                                     )
        callbacks_list = [checkpoint]
        cost = LSTMmodel.fit(train_x, train_y, batch_size=batch_size, epochs=epochs
                         , verbose=1, validation_data=(test_x, test_y), callbacks=callbacks_list)
def tune_ROC_curve(x,y,modelname):
    from sklearn.metrics import roc_auc_score
    from sklearn.metrics import roc_curve
    import tensorflow as tf
    for i in modelname:
        if str(i) =='tbsvm':
            svmx = np.reshape(x, (-1, 200 * 4))
            model = joblib.load(str(i) + '.h5')
            mean_fpr = np.arange(0, 1.1, 0.1)
            score = model.predict_proba(svmx)
            fpr, tpr, thresholds = roc_curve(y, score[:, 1])
            testy_probs = tf.one_hot(y, depth=2)
            lrauc = roc_auc_score(testy_probs, score)
            plt.plot(fpr, tpr,
                     lw=1,
                     alpha=0.3,
                     label='%s(AUC = %0.2f)' % (i, lrauc))

        elif str(i) =='tbLSTM' or str(i)=='tbTCN':
            model = load_model(str(i) + '.h5')
            mean_fpr = np.arange(0, 1.1, 0.1)
            score = model.predict(x)
            fpr, tpr, thresholds = roc_curve(y, score)
            f = (interp1d(fpr, tpr, fill_value="extrapolate"))
            tprs = f(mean_fpr)
            tprs[0] = 0
            lrauc = roc_auc_score(y, score)
            plt.plot(fpr, tpr,
                     lw=1,
                     alpha=0.3,
                     label='%s(AUC = %0.2f)'%(i, lrauc))
        elif str(i) == 'LSTM1' or str(i) == 'TCN1':
            nonx = np.load('x.npy').astype(np.float)
            nony = np.load('y.npy', allow_pickle=True).astype(np.int)
            train_x, test_x, train_y, test_y = train_test_split(nonx, nony, test_size=0.2, random_state=43)
            model = load_model(str(i) + '.h5')
            mean_fpr = np.arange(0, 1.1, 0.1)
            score = model.predict(test_x)
            fpr, tpr, thresholds = roc_curve(test_y, score)
            f = (interp1d(fpr, tpr, fill_value="extrapolate"))
            tprs = f(mean_fpr)
            tprs[0] = 0
            lrauc = roc_auc_score(test_y, score)
            plt.plot(fpr, tpr,
                     lw=1,
                     alpha=0.3,
                     label='%s(AUC = %0.2f)' % (i, lrauc))
        else:
            nonx = np.load('x.npy').astype(np.float)
            nony = np.load('y.npy', allow_pickle=True).astype(np.int)
            train_x, test_x, train_y, test_y = train_test_split(nonx, nony, test_size=0.2, random_state=43)
            svmx = np.reshape(test_x, (-1, 200 * 4))
            model = joblib.load(str(i) + '.h5')
            mean_fpr = np.arange(0, 1.1, 0.1)
            score = model.predict_proba(svmx)
            fpr, tpr, thresholds = roc_curve(test_y, score[:, 1])
            testy_probs = tf.one_hot(test_y, depth=2)
            lrauc = roc_auc_score(testy_probs, score)
            plt.plot(fpr, tpr,
                     lw=1,
                     alpha=0.3,
                     label='%s(AUC = %0.2f)' % (i, lrauc))

    plt.plot([0, 1], [0, 1],
             linestyle='--',
             lw=2,
             color='r',
             label='Chance',
             alpha=.8)
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    # plt.legend(bbox_to_anchor=(1, 0), loc=3, borderaxespad=0)
    plt.legend()
    plt.show()
    plt.clf()
if __name__ == '__main__':
    LSTMss=LSTM_Net()
    LSTMss.train(epochs=50)
    TCNs=TCN_Net()
    TCNs.train(epochs=50)
    train_svm()
    x = np.load('x.npy').astype(np.float)
    y = np.load('y.npy', allow_pickle=True).astype(np.int)
    train_x, test_x, train_y, test_y = train_test_split(x, y, test_size=0.2, random_state=43)
    modelname=['LSTM1','TCN1']
    name2=['tbLSTM','tbTCN','tbsvm']
    tune_ROC_curve(test_x,test_y,modelname)
    tune_ROC_curve(test_x, test_y, name2)

