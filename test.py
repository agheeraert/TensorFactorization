from Tensor import MDTensor
import mdtraj as md
import numpy as np
from tqdm import tqdm
import pickle as pkl
import matplotlib.pyplot as plt
import pandas as pd
from os.path import join as jn
from collections import OrderedDict
from Bio.PDB.Polypeptide import aa1, aa3
import networkx as nx
# import tensorly.decomposition as td
from FactorAnalysis import FactorAnalysis

# def kron_moinsbete(N):
#     M = np.zeros((N,N,N))
#     for i in range(N):
#         M[i,i,i] = 1
#     return M


# Tensor1 = MDTensor('/home/aghee/PDB/prot_apo_sim1_s10.dcd', '/home/aghee/PDB/prot.prmtop')
# Tensor1.load_tensor('results/apo10_sim1.p')
for i in tqdm(range(25,26)):
    Parafac = np.load("results/decal/ABC_decal"+str(i)+".npy", allow_pickle=True)
    A,B,C = Parafac.item()['A'], Parafac.item()['B'], Parafac.item()['C']
    FA = FactorAnalysis(A, B, C)
    FA.plot_membership_distribution('results/decal/membership_distribution_'+str(i)+'.png')
    FA.plot_activity_pattern('results/decal/activity_pattern_'+str(i)+'.png')
    FA.create_labels('/home/aghee/PDB/prot_apo_sim1_s10.dcd', '/home/aghee/PDB/prot.prmtop')
    indexes = []
    for j in range(i):
        indexes.append(FA.single_component_analysis(j,1)[0])
    FA.membership_diagram(indexes, 'results/decal/membership_diagram'+str(i)+'.png')
