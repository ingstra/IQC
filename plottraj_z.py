from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from scipy.misc import factorial
import sympy.mpmath as mp

#import sys
#trace = np.loadtxt('trace.dat')

#traj = np.loadtxt('traj.dat')
#exact = np.loadtxt('exact.dat')
#current = np.loadtxt('current.dat')

#dt = 1e-3
#nruns = 10000

gamma = 1

#mu = 0
#variance = dt*nruns
#sigma = np.sqrt(variance)


#plt.plot(traj[:,0],traj[:,1])
#plt.plot(exact[:,0],exact[:,1],'g--')

#plt.plot(traj[:,0],np.exp(-gamma*traj[:,0]),'r:')

#plt.plot(trace[:,0],trace[:,1])

def p_plus(t):
    Omega = 5
    gamma = 1
    
    if 4*Omega < gamma:
        kappa = np.sqrt(np.power(gamma/4,2)-np.power(Omega,2))
    else:
        kappa = 1j*np.sqrt(-np.power(gamma/4,2)+np.power(Omega,2))
    
    K1 = np.power(Omega,2) / np.add(np.power(gamma,2) , 2*np.power(Omega,2) )
    hyp = np.cosh(kappa*t) + 3*gamma*np.sinh(kappa*t)/(4*kappa)

    return  K1*(1 - np.exp(-3*gamma*t/4)*hyp )

t=np.linspace(0,8,100)

plt.plot(t,p_plus(t),'r--')



plt.tight_layout()

plt.show()
