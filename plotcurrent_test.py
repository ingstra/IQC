from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from scipy.misc import factorial
import sympy.mpmath as mp

current = np.loadtxt('current.dat')
current=current[0,:]
h =  np.loadtxt('test.dat')

dt = 1e-3
nruns = 10000

gamma = 2

mu = 0
variance = dt*nruns
sigma = np.sqrt(variance)

###plt.hist(data, bins=np.arange(min(data), max(data) + binwidth, binwidth))
# values, bins, _ = plt.hist(current,10,normed=0,histtype='bar')


N=50
L=5
dx=2*L/N

calcbin=[-L+n*dx for n in range (0,N)]

plt.hist(current, bins=calcbin)

plt.plot(h[:,0],h[:,1],'*')




plt.show()
