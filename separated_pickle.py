import pickle as pkl
import numpy as np 
import pandas as pd
import os
import re
from rpy2.robjects import r
from rpy2.robjects.numpy2ri import numpy2ri

for i in range(1,5):
    ao = 'results/noH/apo10_sim'+str(i)
    po = 'results/noH/prfar10_sim'+str(i)
    mat = pkl.load(open('results/noH/sim'+str(i)+'.p', 'rb'))
    np.save(ao+'.npy', mat[:,:,:1000:10])
    np.save(po+'.npy', mat[:,:,1000::10])




