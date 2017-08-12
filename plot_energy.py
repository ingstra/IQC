from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from numpy import pi as pi
from scipy.stats import norm
from scipy.misc import factorial

e = np.loadtxt('energy.dat')

plt.plot(e[:,0],e[:,1],'k',linewidth=2,label='a')
print(np.mean(e[:,1]))

#plt.plot(trace[:,0],trace[:,1])

plt.xlabel(r'ns')
plt.ylabel(r'$\omega_{eff}\langle a^\dagger a \rangle$')

plt.tight_layout()
#plt.savefig('testplot',figsize=(20,10))

plt.show()
