from numpy import genfromtxt
import numpy as np
import matplotlib.pyplot as plt

table = genfromtxt('results/noH/partial_cc.csv', delimiter=',')
best5 = np.sort(table, axis=-1)[:,-5:]
plt.style.use('classic')
f = plt.figure()
plt.plot(range(2,best5.shape[0]+2), np.mean(best5, axis=-1), color='k')
plt.plot(range(2,best5.shape[0]+2), [90]*best5.shape[0], color='g', linestyle='--')
for i in range(5):
    plt.scatter(range(2,best5.shape[0]+2), best5[:,i], marker='+', color='b')
plt.xlabel("Number of components")
plt.ylabel("core consistency (%)")
plt.xlim(1.5,best5.shape[0]+1.5)
plt.ylim(0,101)
plt.savefig("results/noH/cc_noH.png")