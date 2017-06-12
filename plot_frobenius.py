from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from scipy.misc import factorial
import sympy.mpmath as mp


h =  np.loadtxt('frobenius.dat')


plt.semilogy(h[:,0],h[:,1])


plt.xlabel('Iteration step')
plt.ylabel('Frobenius norm of matrix difference')

plt.show()
