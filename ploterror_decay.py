from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from scipy.misc import factorial
import sympy.mpmath as mp

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#import sys
#trace = np.loadtxt('trace.dat')

traj = np.loadtxt('traj.dat')
exact = np.loadtxt('exact.dat')
current = np.loadtxt('current.dat')

dt = 1e-3
nruns = 10000

gamma = 1

mu = 0
variance = dt*nruns
sigma = np.sqrt(variance)



#plt.plot(traj[:,0],np.abs(traj[:,1]-exact[:,1]),'g')
plt.plot(traj[:,0],np.abs(traj[:,1]-np.exp(-gamma*traj[:,0])),'k',linewidth=2)
print(np.max(np.abs(traj[:,1]-np.exp(-gamma*traj[:,0]))))
plt.xlabel(r'Time')
plt.ylabel(r'Error')
#plt.plot(exact[:,0],exact[:,1],'g--')

#plt.plot(traj[:,0],np.exp(-gamma*traj[:,0]),'r:')



#plt.plot(trace[:,0],trace[:,1])




def avg_Z(t):
    Omega = 1
    gamma = 1
    
    if 4*Omega < gamma:
        kappa = np.sqrt(np.power(gamma/4,2)-np.power(Omega,2))
    else:
        kappa = 1j*np.sqrt(-np.power(gamma/4,2)+np.power(Omega,2))
    
    K1 = np.power(Omega,2) / np.add(np.power(gamma,2) , 2*np.power(Omega,2) )
    hyp = np.cosh(kappa*t) + 3*gamma*np.sinh(kappa*t)/(4*kappa)

    return  K1*(1 - np.exp(-3*gamma*t/4)*hyp ) #- 1/2

def avg_Z_big_Omega(t):
    Omega = 5
    gamma = 1
    return -np.exp(-3*gamma*t/4)*np.cos(Omega*t)/2

def p_plus(t):
    Omega = 5
    gamma = 1
    return (1 - np.exp(-3*gamma*t/4)*np.cos(Omega*t))/2

#plt.plot(traj[:,0],avg_Z(traj[:,0]),'r--')

#plt.plot(data[:,0],p_plus(data[:,0]),'g')

plt.tight_layout()
plt.savefig('testplot',figsize=(20,10))

plt.show()