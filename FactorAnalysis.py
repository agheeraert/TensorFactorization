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
import seaborn as sns


class FactorAnalysis():
    def __init__(self, A, B, C):
        self.three2one = dict(zip(aa3, aa1))
        self.N, self.R = [*A.shape]
        self.T = C.shape[0]
        self.A = A
        self.C = C
    
    def create_labels(self, trajs, topo):
        traj = md.load(trajs, top=topo)
        topo = traj.topology.to_dataframe()[0]
        resId2resName = pd.Series(topo['resName'].values, topo['resSeq'].values).to_dict()
        self.resId2resName = {}
        for elt in resId2resName:
            if elt <= 253:
                self.resId2resName[elt] = self.three2one[resId2resName[elt]]+str(elt+1)+':F'
            else:
                self.resId2resName[elt] = self.three2one[resId2resName[elt]]+str(elt-253)+':H'
    
    def save_labels(self, output):
        if self.resId2resName:
            pkl.dump(self.resId2resName, open(output, "wb"))
        else:
            print("No labels created")
    
    def load_labels(self, path):
        self.resId2resName = pkl.load(open(path, "rb"))

    def plot_membership_distribution(self, output):
        f = plt.figure()
        for component in range(self.R):
            plt.subplot(self.R//3+1-(self.R%3==0),3,component+1)
            plt.title('Component '+str(component+1))
            plt.hist(self.A[:,component], color='k', log=True)
            plt.xlabel('Level of membership')
            if component % 3 == 0:
                plt.ylabel('#occurences)')
        plt.tight_layout()
        plt.savefig(output)
        plt.close()
    
    def activity_pattern(self):
        return np.multiply(self.C, np.sum(self.A, axis=0))
    
    def plot_activity_pattern(self, output):
        S = self.activity_pattern()
        f= plt.figure()
        for component in range(self.R):
            plt.subplot(self.R//3+1-(self.R%3==0),3,component+1)
            plt.title('Component '+str(component+1))
            plt.plot(S[:,component], color='k')
            plt.xlabel('Time (ns)')
            if component % 3 == 0:
                plt.ylabel('Strength')
        plt.tight_layout()
        plt.savefig(output)
        plt.close()                    
