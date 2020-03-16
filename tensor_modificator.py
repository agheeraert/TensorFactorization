import numpy as np
from tqdm import tqdm
import pickle as pkl
import matplotlib.pyplot as plt
import pandas as pd
from os.path import join as jn
from scipy.spatial.distance import pdist, squareform
from collections import OrderedDict
from Bio.PDB.Polypeptide import aa1, aa3
import networkx as nx
import argparse
import mdtraj as md
from Tensor import Tensor

parser = argparse.ArgumentParser(description='Operates on ')
parser.add_argument('-f',  type=str, nargs='+',
                    help='List of path of tensors to concatenate and modify')
parser.add_argument('-w', type=int, nargs=1, default=10, 
                    help='time window for perfoming the operations')
parser.add_argument('-a',  type=str, default='all',
                    help='actions to perform')
parser.add_argument('-o',  type=str, default='results/diff.png',
                    help='output path')
args = parser.parse_args()

ts = Tensor(args.f)
ts.set_window(args.w)
if args.a == 'all':
    ts.all_operations(args.o)



    

