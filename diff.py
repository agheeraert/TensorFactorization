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

parser = argparse.ArgumentParser(description='Plot cc from R npy output(s)')
parser.add_argument('f',  type=str, nargs='+',
                    help='list of trajectories and topology file')
parser.add_argument('-o',  type=str, default='results/diff.png',
                    help='output path')
args = parser.parse_args()

topo = args.f[-1]
trajs = args.f[:-1]
diff = []
traj = md.load(trajs, top=topo)
for i in tqdm(range(traj.n_frames-1)):
    diff.append(md.rmsd(traj[i], traj[i+1]))

f = plt.figure()        
plt.plot(range(1, traj.n_frames), diff, color='k')
plt.title('Diff vs frame number')
plt.xlabel('Frame number')
plt.ylabel('rmsd')
plt.savefig(args.o)
