from numpy import genfromtxt
import numpy as np
import matplotlib.pyplot as plt

table = genfromtxt('results/partial_cc.csv', delimiter=',')
best5 = np.sort(table, axis=-1)[:,-5:]

f = plt.figure()
plt.plot(range(2,13), np.mean(best5, axis=-1), color='k')
for i in range(5):
    plt.scatter(range(2,13), best5[:,i], marker='+', color='g')
plt.xlabel("Number of components")
plt.ylabel("core consistency (%)")
plt.ylim(0,100)
plt.savefig("results/coreconsistency.png")