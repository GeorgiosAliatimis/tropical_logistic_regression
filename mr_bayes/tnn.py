import tensorflow as tf
from tensorflow.keras.layers import Layer, Dense
from tensorflow.keras.models import Sequential
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.backend import repeat_elements
import numpy as np
from sklearn.metrics import auc, roc_curve
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys

class TropEmbed(Layer):
    def __init__(self, units=2, input_dim=3):
        super(TropEmbed, self).__init__()
        self.w = self.add_weight(shape=(units, input_dim), \
                                 initializer="random_normal", trainable=True)
        self.units = units
        self.input_dim = input_dim
    def call(self, x):
        x_reshaped = tf.reshape(x,[-1, 1, self.input_dim])
        x_for_broadcast = repeat_elements(x_reshaped, self.units, 1)
        values, indices = tf.math.top_k(x_for_broadcast + self.w, self.input_dim)
        return values[:,:,0] - values[:,:,-1] # symmetric tropical distance
      
class Normalize(Layer):
    def __init__(self):
        super(Normalize,self).__init__()
        self.w = self.add_weight(shape=(1,), initializer="ones", trainable=True)
    def call(self,inputs):
        return 10*self.w * inputs

def train_model(x_train,y_train,hidden_layer):
    #model = Sequential([Normalize(),hidden_layer,Dense(1, activation="sigmoid")])
    model = Sequential([hidden_layer,Dense(1, activation="sigmoid")])
    model.compile(optimizer=Adam(0.1), loss="binary_crossentropy", metrics=["accuracy"])
    model.fit(x_train, y_train, epochs=20)
    return model

def train_tropical_model(x_train,y_train,num_neurons=2):
    tree_dim = x_train.shape[1] 
    hidden_layer = TropEmbed(num_neurons, tree_dim)
    return train_model(x_train,y_train,hidden_layer)

def get_accuracy(y_test_pred,y_test):
    p_thres = 0.5
    acc = ((y_test_pred > p_thres).reshape([-1,]) == y_test.reshape([-1,])).mean()
    return acc 

def get_auc(y_test_pred,y_test,plot_roc=False):
    fpr, tpr, thresholds = roc_curve(y_test, y_test_pred)
    if plot_roc:    
        plt.plot(fpr,tpr)
        plt.show()
    return auc(fpr, tpr)

def classify_to_methods(D1,D2,model="tropical"):
    D1 = D1[0::2,:]
    D2 = D2[0::2,:]
    len1 = D1.shape[0]
    len2 = D2.shape[0]
    x = np.vstack([D1,D2])
    y = np.hstack([np.repeat(0, len1), np.repeat(1, len2)]).reshape((len1+len2,1))
    l = np.arange(0, len1+len2)
    np.random.shuffle(l)
    train_size = int(0.8 * (len1+len2))
    train = l[0:train_size]
    test = l[train_size:]
    x_train = x[train,]
    x_test = x[test,]
    y_train = y[train]
    y_test = y[test]
    if model == "tropical":
        trained_model = train_tropical_model(x_train,y_train)
    elif model == "classical":
        trained_model = train_classical_model(x_train,y_train)
    else:
        raise Exception("Method not found")
    y_test_pred = trained_model.predict(x_test)
    aucval = get_auc(y_test_pred,y_test)
    print(f"AUC={aucval}")
    #print(f"norm={trained_model.layers[0].get_weights()[0]}")
    return aucval

files = sorted(os.listdir("trees1"), key=lambda x: int(x.split(".")[0]))
aucs=list()
l = float(sys.argv[1])
for file in files:
    D_1 =  l*np.loadtxt(open(f"trees1/{file}","rb"),delimiter=" ")
    D_2 =  l*np.loadtxt(open(f"trees2/{file}","rb"),delimiter=" ")
    s = [classify_to_methods(D_1,D_2,model="tropical") for _ in range(5)]
    aucval = float(np.median(s))
    print(file, aucval)
    aucs.append(aucval)
aucs = np.array(aucs)

print(aucs)
# plt.plot(aucs)
# plt.show()

with open('aucs','w') as f:
	f.write('\n'.join(map(str,aucs)))