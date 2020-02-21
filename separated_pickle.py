import pickle as pkl
import numpy as np 
import pandas as pd
import os
import re
from rpy2.robjects import r
from rpy2.robjects.numpy2ri import numpy2ri

for i in range(1,5):
    ao = 'results/apo10_sim'+str(i)
    po = 'results/prfar10_sim'+str(i)
    mat = pkl.load(open('results/sim'+str(i)+'.p', 'rb'))
    pkl.dump(mat[:,:,:1000:10], open(ao+'.p', 'wb'))
    pkl.dump(mat[:,:,1000::10], open(po+'.p', 'wb'))




