import numpy as np
import matplotlib.pyplot as plt

# read data
processes = [1,2,4]
time = []
for pro in processes:
    time.append(int(np.loadtxt('time'+str(pro)+'.dat')))
plt.plot(processes, time, 'bo')
plt.xlabel('# of Processes')
plt.ylabel('Time (seg)')

plt.savefig('time.png')
