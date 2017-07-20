from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from numpy import pi as pi
from scipy.stats import norm
from scipy.misc import factorial



traj = np.loadtxt('noise.dat')

#plt.plot(traj[:,0],traj[:,1],'b',linewidth=2)

values, bins, _ = plt.hist(traj[:,1],50)

mufit,stdfit = norm.fit(traj[:,1])

print 'mean: ',mufit
print 'std: ', stdfit



plt.tight_layout()
#plt.savefig('testplot',figsize=(20,10))

plt.show()
