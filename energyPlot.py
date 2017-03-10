import numpy as np
import matplotlib.pyplot as plt

P = 1000
k = 3
data = np.loadtxt('energy.dat')
data = data.reshape(k, P)

plt.plot(data.T)
plt.savefig('energy.png')
