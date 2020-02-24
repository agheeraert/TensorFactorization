from os.path import join as jn
import numpy as np
import argparse
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='Plot cc from R npy output(s)')
parser.add_argument('f',  type=str, nargs='+',
                    help='Numpy file containing the core consistencies')
args = parser.parse_args()

for fichier in args.f:
    table = np.load(fichier)
    if '.npy' in fichier:
        output = fichier.replace('.npy', '.png')
    else:
        output = fichier+'.png'
    best5 = np.sort(table, axis=-1)[:,-5:]
    f = plt.figure()
    plt.plot(range(2,best5.shape[0]+2), np.mean(best5, axis=-1), color='k')
    for i in range(5):
        plt.scatter(range(2,best5.shape[0]+2), best5[:,i], marker='+', color='g')
    plt.xlabel("Number of components")
    plt.ylabel("core consistency (%)")
    plt.ylim(0,100)
    plt.savefig(output)

