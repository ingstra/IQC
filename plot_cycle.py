from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from numpy import pi as pi
from scipy.stats import norm
from scipy.misc import factorial

cycle = np.loadtxt('cycle.dat')

plt.plot(cycle[:,0],np.log(cycle[:,1]),'k',linewidth=2)



#plt.plot(trace[:,0],trace[:,1])

plt.xlabel(r'$\omega_{eff}/\omega_a$')
plt.ylabel(r'$\log(U_a/\omega_a)$')

plt.tight_layout()
#plt.savefig('testplot',figsize=(20,10))

plt.show()
