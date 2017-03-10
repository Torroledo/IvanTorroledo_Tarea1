import numpy as np
import matplotlib.pyplot as plt

# read data
processors = [1,2,4]
time = []
for pro in processors:
    time.append(int(np.loadtxt('time'+str(pro)+'.dat')))
plt.plot(processors, time, 'ro')
plt.xlabel('# of Processors')
plt.ylabel('Time ($seg$)')
plt.savefig('time.png')
