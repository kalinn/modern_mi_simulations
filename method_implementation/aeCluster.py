#Python 3.8.1 

import os
import pandas as pd
import numpy as np
import keras
from keras import backend as K
from sklearn import preprocessing
import sys
sys.path.append("/project/flatiron_ucc")
sys.path.append("/project/flatiron_ucc/__pycache__")
import desc
import tensorflow
import session_info
os.chdir("/project/flatiron_ucc")
rootdir = "/project/flatiron_ucc/programs/kylie/RunMe2/final_results"
#if os.path.isdir(rootdir)==False:
#    os.mkdir (rootdir)

mechList = ("MCAR","MAR","MNAR1","MNAR2")
proportionList = (10,30,50)
for w in range(0,len(mechList)):
    mech = mechList[w]
    foldName = "AE" + mech
    path = rootdir + "/" + foldName
#    if os.path.isdir (path)==False:
#        os.mkdir (path)
    for k in range(0,len(proportionList)):
        propName = proportionList[k]
        path = rootdir + "/" + foldName + "/" + str (propName)
#        if os.path.isdir (path)==False:
#            os.mkdir (path)
        i = sys.argv[1]
        filename = "".join(("programs/kylie/RunMe2/datasets/mDats/",mech,"/",str(propName),"/","mData",str(i),".csv"))
        train = pd.read_csv(filename)
        time = train["time"]
        train = train.drop("time",axis=1)
        #16 nodes in 2nd layer
        n = len(train.columns)
        dim = [train.shape[1],n+7,n+14,n+21]
        #dim= [train.shape[1],n+7,n+14,n+21,n+14,n+7]
        #dim= [train.shape[1],n-5,n-10]
        trainFilled = train.fillna(train.mean()) #fill NAs with column mean b/c encoder needs full dataset
        train2 = np.matrix(trainFilled)
        for x in range(1, 11):
            K.clear_session()
            try:
                sae = desc.SAE(dims=dim,drop_rate=0.5,batch_size=256,random_seed=x,act="tanh",actincenter="linear")
            except:
                print ("desc.SAE error")
            finally:
                sae = desc.SAE(dims=dim,drop_rate=0.5,batch_size=256,random_seed=x,act="tanh",actincenter="linear")
            sae.fit(train2,epochs=500)
            predict = sae.autoencoders.predict(train2)
            #predict1= predict.round(0)
            predictNew = pd.DataFrame(predict,columns=train.columns,index=train.index)
            miss = pd.DataFrame(train.copy())
            miss[miss.isnull()] = predictNew
            final = miss
            final.insert(loc=2,column="time",value=time)
            final = pd.DataFrame(final)
            filename = "".join(("programs/kylie/RunMe2/final_results/","AE",mech,"/",str(propName),"/","result",str(i),"_num_",str(x),".csv"))
            final.to_csv(filename)


