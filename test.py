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
from Tensor import Tensor

# def kron_moinsbete(N):
#     M = np.zeros((N,N,N))
#     for i in range(N):
#         M[i,i,i] = 1
#     return M


# Tensor1 = MDTensor('/home/aghee/PDB/prot_apo_sim1_s10.dcd', '/home/aghee/PDB/prot.prmtop')
# Tensor1.load_tensor('results/apo10_sim1.p')


for i in tqdm(range(8,20)):
    Parafac = np.load("results/shuffled_rolling/ABC"+str(i)+".npy", allow_pickle=True)
    A,B,C = Parafac.item()['A'], Parafac.item()['B'], Parafac.item()['C']
    FA = FactorAnalysis(A, B, C)
    FA.reorder('results/shuffled_rolling/apo_1_index.npy')
    FA.plot_membership_distribution('results/shuffled_rolling/membership_distribution_'+str(i)+'.png')
    FA.plot_activity_pattern('results/shuffled_rolling/activity_pattern_'+str(i)+'.png')
    FA.create_labels('/home/aghee/PDB/prot_apo_sim1_s10.dcd', '/home/aghee/PDB/prot.prmtop')
    indexes = []
    for j in range(i):
        indexes.append(FA.single_component_analysis(j,2.5)[0])
    FA.membership_diagram(indexes, 'results/shuffled_rolling/membership_diagram'+str(i)+'.png')
    FA.components_to_vmd(indexes, 'results/shuffled_rolling/'+str(i)+'.tcl')

# tensor = Tensor(['results/apo_sim1.npy', 'results/prfar_sim1.npy'])
# tensor.set_window(10)
# tensor.all_operations('results/sim1/')
