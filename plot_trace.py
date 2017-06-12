from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from scipy.misc import factorial
import sympy.mpmath as mp

#import sys
trace = np.loadtxt('trace.dat')


plt.plot(trace[:,0],trace[:,1])



plt.show()
