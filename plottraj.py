from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from numpy import pi as pi
from scipy.stats import norm
from scipy.misc import factorial

#import sys
#trace = np.loadtxt('trace.dat')

traj = np.loadtxt('traj.dat')
#exact = np.loadtxt('exact.dat')
#current = np.loadtxt('current.dat')


nh=0.2
wa = 2*pi*10
wb = 2*pi*0.5
ka=2*pi*2
kb = 2*pi*0.05
g=2*pi*0.5

#dnasquared=nh*(nh+1)
#C=4*g**2/(ka*kb)
#Cprime = C/(1+4*(wb/(4*ka))**2)
#Pss= np.ones_like(traj[:,0])*kb*wb*Cprime*dnasquared/2

#print(Cprime)
#print(Pss[1])
#plt.plot(traj[:,0],Pss,'r--',linewidth=2)

plt.plot(traj[:,0],traj[:,1],'k',linewidth=2)
#plt.plot(exact[:,0],exact[:,1],'g--',linewidth=2)


#plt.plot(trace[:,0],trace[:,1])

plt.xlabel(r'$\omega_a\,t$')
plt.ylabel(r'$\langle n_b \rangle$')

plt.tight_layout()
#plt.savefig('testplot',figsize=(20,10))

plt.show()
