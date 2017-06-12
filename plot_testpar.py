from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from scipy.misc import factorial
import sympy.mpmath as mp

current = np.loadtxt('current.dat')

nangles = len(current)


hist=[ np.histogram(current[i],50,normed=1) for i in range(0,nangles)]

values = [hist[i][0] for i in range(0,nangles)]


bins = [hist[i][1] for i in range(0,nangles)]


bincenters = [ 0.5*(bins[i][1:]+bins[i][:-1]) for i in range(0,nangles)]


plt.plot(bincenters,values,'ko')

range = 5
x= np.linspace(-range,range)
h_poly = np.frompyfunc(mp.hermite,2,1)
n=0
const0 = 1/(np.sqrt(2**n * factorial(n)) * (2*np.pi)**0.25 )
rest0 = const0* np.exp(-x**2/4)
tot0  = (rest0*h_poly(n,x/np.sqrt(2)))**2
plt.plot(x,tot0,'r-')

plt.show()
