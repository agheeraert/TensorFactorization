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

parser = argparse.ArgumentParser(description='Operates on ')
parser.add_argument('f',  type=str, nargs='+',
                    help='list of trajectories and topology file')
parser.add_argument('-p',  type=str, 
                    help='output path')
parser.add_argument('-o',  type=str, default='results/diff.png',
                    help='output path')
args = parser.parse_args()

if 'all' in args.p:
    

