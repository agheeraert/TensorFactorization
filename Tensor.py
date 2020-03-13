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
# import tensorly.decomposition as td
# from tensorly.tucker_tensor import tucker_to_tensor
# from sktensor import dtensor, cp_als
# from tensorly import to_numpy as tn

class MDTensor():
    def __init__(self, trajs, topo, cutoff=5, discard_hydrogens=False):
        """
        cutoff in ANGSTROM"""
        self.traj = md.load(trajs, top=topo)
        if discard_hydrogens:
            self.traj = self.discard_hydrogens()
        print(self.traj)
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
        topo = self.traj.topology.to_dataframe()[0]
        isH = np.char.startswith(np.array(topo['name'], dtype=str), 'H')
        indices = np.where(isH==0)[0]
        return self.traj.restrict_atoms(indices)

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

class Tensor():
    def __init__(self, paths):
        self.tensor = np.concatenate([np.load(path) for path in paths], axis=-1)
        self.n_nodes = self.tensor.shape[0]
        self.n_frames = self.tensor.shape[-1]

    def set_window(self, window):
        self.window = window

    def shuffle(self, output):
        perm = np.random.permutation(range(self.n_nodes))
        self.tensor = self.tensor[perm]
        self.tensor = self.tensor[:,perm]
        np.save(output, self.tensor)
        assert '.npy' in output, 'Output must be in .npy format' 
        np.save(output.replace('.npy', '_index.npy'), perm)

    def _mean(self, tensor):
        return np.round(np.mean(tensor.reshape(self.n_nodes, self.n_nodes, -1, self.window), axis=-1))
    
    def mean(self, output):
        np.save(output, self._mean(self.tensor))
    
    def decmean(self, output):
        mean = self._mean(self.tensor)
        mean5 = self._mean(np.roll(self.tensor, self.window//2))
        decmean = np.zeros((self.n_nodes, self.n_nodes, self.n_frames//self.window*2))
        decmean[:,:,::2] = mean
        decmean[:,:,1::2] = mean5
        np.save(output, decmean)
    
    def rollmean(self, output):
        roll = np.zeros((self.n_nodes, self.n_nodes, self.n_frames))
        arrays = [self._mean(np.roll(self.tensor, i)) for i in range(self.window)]
        for j, elt in enumerate(arrays):
            for i in range(self.n_frames//self.window):
                roll[:,:,i*10+j] = elt[:,:,i]
        np.save(output, roll)
    
    def all_operations(self, folder):
        self.shuffle(jn(folder, 'shuffle.npy'))
        self.mean(jn(folder, 'mean.npy'))
        self.decmean(jn(folder, 'decmean.npy'))
        self.rollmean(jn(folder, 'rollmean.npy'))


if __name__ == '__main__':
    for i  in range(1,2):
        Tensor1 = MDTensor(['/home/aghee/PDB/prot_apo_sim'+str(i)+'_s10.dcd', '/home/aghee/PDB/prot_prfar_sim'+str(i)+'_s10.dcd'], '/home/aghee/PDB/prot.prmtop', discard_hydrogens=True)
        Tensor1.create_tensor()
        Tensor1.save('results/noH/sim'+str(i)+'.p')
        del Tensor1
    # Tensor1 = MDTensor(['/home/aghee/PDB/prot_apo_sim1_s10.dcd', '/home/aghee/PDB/prot_prfar_sim1_s10.dcd'], '/home/aghee/PDB/prot.prmtop')
    # Tensor1.load_tensor('results/sim1.p')
    # Tensor1.assess_quality()
