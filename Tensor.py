import mdtraj as md
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
import tensorly.decomposition as td
from tensorly.tucker_tensor import tucker_to_tensor
from sktensor import dtensor, cp_als
from tensorly import to_numpy as tn

class MDTensor():
    def __init__(self, trajs, topo, cutoff=5, discard_hydrogens=False):
        """
        cutoff in ANGSTROM"""
        self.traj = md.load(trajs, top=topo)
        if discard_hydrogens:
            self.traj = self.discard_hydrogens()
        self.cutoff = cutoff/10 #conversion in nm
        topo = self.traj.topology.to_dataframe()[0]
        self.atom2res = pd.Series(topo['resSeq'].values, topo.index).to_dict()
        self.three2one = dict(zip(aa3, aa1))
        resId2resName = pd.Series(topo['resName'].values, topo['resSeq'].values).to_dict()
        self.resId2resName = {}
        for elt in resId2resName:
            if elt <= 253:
                self.resId2resName[elt] = self.three2one[resId2resName[elt]]+str(elt+1)+':F'
            else:
                self.resId2resName[elt] = self.three2one[resId2resName[elt]]+str(elt-253)+':H'
    
    def discard_hydrogens(self):
        return self.traj #TO DO

    def save(self, output):
        pkl.dump(self.tensor, open(output, 'wb'))
    
    def create_tensor(self):
        L_tensors = []
        for i in tqdm(range(self.traj.n_frames)):
            distances = pdist(self.traj.xyz[i])
            indexes = np.where(squareform(distances)<=self.cutoff)
            u, v = np.vectorize(self.atom2res.get)(indexes)
            net = nx.Graph()
            net.add_edges_from((u[i], v[i]) for i in range(len(u)))
            net = nx.relabel_nodes(net, self.resId2resName)
            L_tensors.append(nx.to_numpy_array(net, dtype=bool))
        self.tensor = np.stack(L_tensors, axis=-1)
    
    def load_tensor(self, file):
        self.tensor = pkl.load(open(file, "rb"))
    
    def factorize(self, R):
        return [*td.non_negative_parafac(self.tensor.astype("float64"), R)]



if __name__ == '__main__':
    for i  in range(2,5):
        Tensor1 = MDTensor(['/home/aghee/PDB/prot_apo_sim'+str(i)+'_s10.dcd', '/home/aghee/PDB/prot_prfar_sim'+str(i)+'_s10.dcd'], '/home/aghee/PDB/prot.prmtop')
        Tensor1.create_tensor()
        Tensor1.save('results/sim'+str(i)+'.p')
        del Tensor1
    # Tensor1 = MDTensor(['/home/aghee/PDB/prot_apo_sim1_s10.dcd', '/home/aghee/PDB/prot_prfar_sim1_s10.dcd'], '/home/aghee/PDB/prot.prmtop')
    # Tensor1.load_tensor('results/sim1.p')
    # Tensor1.assess_quality()