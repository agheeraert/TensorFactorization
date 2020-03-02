from numpy import genfromtxt
import numpy as np
import matplotlib.pyplot as plt

table = genfromtxt('results/partial_cc.csv', delimiter=',')
best5 = np.sort(table, axis=-1)[:,-5:]
plt.style.use('classic')
f = plt.figure()
plt.plot(range(2,best5.shape[0]+2), np.mean(best5, axis=-1), color='k')
for i in range(5):
    plt.scatter(range(2,best5.shape[0]+2), best5[:,i], marker='+', color='g')
plt.xlabel("Number of components")
plt.ylabel("core consistency (%)")
plt.xlim(1.5,best5.shape[0]+1.5)
plt.ylim(0,100)
plt.savefig("results/cc_classic.png")