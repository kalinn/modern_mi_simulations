import sys
sys.path.append("/project/flatiron_ucc/programs/kylie/RunMe4/code_kal")

import os
import random
import numpy as np
import pandas as pd
from math import sqrt

from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import MinMaxScaler

import torch
import torch.nn as nn
import torch.utils.data
import torch.optim as optim

os.chdir("/project/flatiron_ucc")

j = sys.argv[1]

mech = str(sys.argv[2])
foldName = "AE" + mech

propName = sys.argv[3]

rootdir = "/project/flatiron_ucc/programs/kylie/RunMe4/final_lowDim"
fp = "programs/kylie/RunMe4/"

data_path = "".join((fp, "datasets/mDats/", mech, "/", str(propName), "/", "mData", str(j), ".csv"))

class Autoencoder(nn.Module):
    def __init__(self, dim):
        super(Autoencoder, self).__init__()
        self.dim = dim
        
        self.drop_out = nn.Dropout(p=0.5)
        
        self.encoder = nn.Sequential(
            nn.Linear(dim+theta*0, dim+theta*1),
            nn.Tanh(),
            nn.Linear(dim+theta*1, dim+theta*2)
        )
            
        self.decoder = nn.Sequential(
            nn.Linear(dim+theta*2, dim+theta*1),
            nn.Tanh(),
            nn.Linear(dim+theta*1, dim+theta*0)
        )
        
    def forward(self, x):
        x = x.view(-1, self.dim)
        x_missed = self.drop_out(x)
        
        z = self.encoder(x_missed)
        out = self.decoder(z)
        
        out = out.view(-1, self.dim)
        
        return out

# Tuning params
theta = 4 # How many nodes to add at each layer of encoder
num_epochs = 500
dropout_ratio = 0.5

use_cuda = False
batch_size  = 256 # not in the paper

dat = pd.read_csv(data_path)
time = dat["time"]
dat = dat.drop("time", axis=1)

complete = dat.dropna()
rows, cols = complete.shape
train_data = complete.values

# standardized between 0 and 1
scaler = MinMaxScaler()
scaler.fit(train_data)
train_data = scaler.transform(train_data)

# Set missing to 0
normd = scaler.transform(dat.values)
normd = pd.DataFrame(normd)
nomiss = normd.fillna(0)
filled_data = nomiss.copy()
filled_data = filled_data.values
filled_data = torch.from_numpy(filled_data).float()

train_data = torch.from_numpy(train_data).float()
train_loader = torch.utils.data.DataLoader(dataset=train_data,batch_size=batch_size,shuffle=True)

# Define Model
device = torch.device("cuda" if use_cuda else "cpu")
model = Autoencoder(dim=cols).to(device)

# Define Loss and Optimizer
loss = nn.MSELoss()
optimizer = optim.SGD(model.parameters(), momentum=0.99, lr=0.01, nesterov=True)

# Train Model
cost_list = []
early_stop = True

for x in range(1, 11):
    torch.manual_seed(x)
    random.seed(x)

    for epoch in range(num_epochs):
        
        total_batch = len(train_data) // batch_size
        
        for i, batch_data in enumerate(train_loader):
            
            batch_data = batch_data.to(device)
            
            reconst_data = model(batch_data)
            cost = loss(reconst_data, batch_data)
            
            optimizer.zero_grad()
            cost.backward()
            optimizer.step()
                    
            if (i+1) % (total_batch//2) == 0:
                print('Epoch [%d/%d], lter [%d/%d], Loss: %.6f'
                    %(epoch+1, num_epochs, i+1, total_batch, cost.item()))
                
            # early stopping rule 1 : MSE < 1e-06
            if cost.item() < 1e-06 :
                early_stop = True
                break

            cost_list.append(cost.item())

        if early_stop :
            break
            
    print("Learning Finished!")

    # Test Model
    model.eval()
    imputed_data = model(filled_data.to(device))
    imputed_data = imputed_data.cpu().detach().numpy()
    final_imp = scaler.inverse_transform(imputed_data)

    final = pd.DataFrame(dat.copy())
    final[final.isnull()] = final_imp
    final.insert(loc=2, column="time", value=time)
    final = pd.DataFrame(final)
    filename = "".join((fp, "final_lowDim/", "AE", mech, "/", str(propName), "/", "result", str(j), "_num_", str(x), ".csv"))
    final.to_csv(filename)


