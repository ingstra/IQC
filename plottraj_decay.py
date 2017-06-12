from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from scipy.misc import factorial

#import sys
#trace = np.loadtxt('trace.dat')

traj = np.loadtxt('traj.dat')
#exact = np.loadtxt('exact.dat')
#current = np.loadtxt('current.dat')



plt.plot(traj[:,0],traj[:,1],'k',linewidth=2)
#plt.plot(exact[:,0],exact[:,1],'g--',linewidth=2)


#plt.plot(trace[:,0],trace[:,1])

plt.xlabel(r'$\omega_a\,t$')
plt.ylabel(r'$\langle n_b \rangle$')

plt.tight_layout()
#plt.savefig('testplot',figsize=(20,10))

plt.show()
