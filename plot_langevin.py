from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from numpy import pi as pi
from scipy.stats import norm
from scipy.misc import factorial

langevin = np.loadtxt('langevin.dat')

#plt.plot(langevin[:,0],langevin[:,1],'k',linewidth=2,label='a')
plt.plot(langevin[:,0],langevin[:,2],'r',linewidth=2,label='b')

plt.legend()

#plt.plot(trace[:,0],trace[:,1])

plt.xlabel(r'ns')
#plt.ylabel(r'$\log(U_a/\omega_a)$')

plt.tight_layout()
#plt.savefig('testplot',figsize=(20,10))

plt.show()
