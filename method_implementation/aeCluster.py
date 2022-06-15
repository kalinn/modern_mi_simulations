#Python 3.8.1 

import os
import pandas as pd
import numpy as np
import keras
from keras import backend as K
from sklearn import preprocessing
import sys
sys.path.append("/project/flatiron_ucc/programs/kylie/RunMe4/code_kal")
sys.path.append("/project/flatiron_ucc/programs/kylie/RunMe4/code_kal/__pycache__")
import desc
import tensorflow
import session_info
os.chdir("/project/flatiron_ucc")
simple = str(sys.argv[3])
if simple=='TRUE':
    rootdir = "/project/flatiron_ucc/programs/kylie/RunMe4/final_simple"
else:
    rootdir = "/project/flatiron_ucc/programs/kylie/RunMe4/final_results"
fp = "programs/kylie/RunMe4/"
#if os.path.isdir(rootdir)==False:
#    os.mkdir (rootdir)

mech = str(sys.argv[2])
foldName = "AE" + mech
path = rootdir + "/" + foldName

proportionList = (10, 30, 50)
for k in range(0, len(proportionList)):
    propName = proportionList[k]
    print(propName)
    path = rootdir + "/" + foldName + "/" + str(propName)
    i = sys.argv[1]
    print(i)
    filename = "".join((fp, "datasets/mDats/", mech, "/", str(propName), "/", "mData", str(i), ".csv"))
    train = pd.read_csv(filename)
    time = train["time"]
    train = train.drop("time", axis=1)
    trainComplete = train.dropna(axis=0)
    #16 nodes in 2nd layer
    n = len(trainComplete.columns)
    dim = [trainComplete.shape[1], n+7, n+14, n+21]
    #dim= [train.shape[1],n-5,n-10]
    # fill NAs with column mean b/c encoder needs full dataset
    trainFilled = np.matrix (train.fillna(value=0))
    train2 = np.matrix(trainComplete)
    for x in range(1, 11):
#        np.random.RandomState(seed=x)
#        N = train2.shape[0]
#        p = train2.shape[1]
#        R = round (N*p*0.1)
#        arr = np.array([None] * R + [1] * (N*p-R))
#        np.random.shuffle(arr)
#        oparr = pd.DataFrame(arr.reshape((N, p)))
#        train2[oparr.isnull()] = 0
        K.clear_session()
        try:
            sae = desc.SAE(dims=dim, drop_rate=0.2, batch_size=32, random_seed=x, act="tanh", actincenter="linear")
        except:
            print("desc.SAE error")
        finally:
            sae = desc.SAE(dims=dim, drop_rate=0.5, batch_size=256, random_seed=x, act="tanh", actincenter="linear")
        sae.fit(train2, epochs=500)
#        predict = sae.autoencoders.predict(trainFilled)
        predict = sae.extract_preds(trainFilled)
        #predict1= predict.round(0)
        predictNew = pd.DataFrame(predict, columns=trainFilled.columns, index=trainFilled.index)
        miss = pd.DataFrame(train.copy())
        miss[miss.isnull()] = predictNew
        final = miss
        final.insert(loc=2, column="time", value=time)
        final = pd.DataFrame(final)
        if simple=='TRUE':
            filename = "".join((fp, "final_simple/", "AE", mech, "/", str(propName), "/", "result", str(i), "_num_", str(x), ".csv"))
        else:
            filename = "".join((fp, "final_results/", "AE", mech, "/", str(propName), "/", "result", str(i), "_num_", str(x), ".csv"))
        final.to_csv(filename)
