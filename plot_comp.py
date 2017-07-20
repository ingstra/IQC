from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from numpy import pi as pi
from scipy.stats import norm
from scipy.misc import factorial



traj = np.loadtxt('traj.dat')

plt.plot(traj[:,0],traj[:,2],'b',linewidth=2)
plt.plot(traj[:,0],traj[:,1],'b',linewidth=2,label=r'quantum')


langevin = np.loadtxt('langevin.dat')

plt.plot(langevin[:,0],langevin[:,1],'--k',linewidth=2)
plt.plot(langevin[:,0],langevin[:,2],'--k',linewidth=2,label='classical')

plt.legend()

plt.xlabel(r'ns')


plt.tight_layout()
#plt.savefig('testplot',figsize=(20,10))

plt.show()
