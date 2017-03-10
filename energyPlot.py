import numpy as np
import matplotlib.pyplot as plt

P = 1000
k = 3
N = 64
delta_t = 5E-3
T = 5.0*(N**2.2)
t = np.linspace(0.0000,T, P)

data = np.loadtxt('energy.dat')
data = data.reshape(k, P)

fig = plt.figure()
ax = plt.subplot(111)
for i in xrange(k):
    number = i+1
    ax.plot(t, data[i,:].T, label='$k = %i$' % number)

plt.xlabel('Time')
plt.ylabel('Energy')
ax.legend()
plt.savefig('energy.png')
