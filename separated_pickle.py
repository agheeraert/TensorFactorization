import pickle as pkl
import numpy as np 
import pandas as pd
import os
import re
from rpy2.robjects import r
from rpy2.robjects.numpy2ri import numpy2ri

interval = 10

for i in range(1,5):
    ao = 'results/noH/apo_mean'+str(interval)+'_sim'+str(i)
    po = 'results/noH/prfar_mean'+str(interval)+'_sim'+str(i)
    mat = pkl.load(open('results/noH/sim'+str(i)+'.p', 'rb'))
    newmat = np.mean(mat.reshape((*mat.shape[:2], int(2000/interval), interval)), axis=-1)
    np.save(ao+'.npy', newmat[:,:,:100])
    np.save(po+'.npy', newmat[:,:,100:])
    




