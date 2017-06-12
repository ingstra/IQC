from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from scipy.misc import factorial
import sympy.mpmath as mp

current = np.loadtxt('current.dat')

dt = 1e-3
nruns = 1300

gamma = 2

mu = 0
variance = dt*nruns
sigma = np.sqrt(variance)

values, bins, _ = plt.hist(current,50,normed=1,histtype='bar')
area = sum(np.diff(bins)*values)
#print bins

mufit, stdfit = norm.fit(current)
#print 'fitted mean: ',mufit #, 'calculated mean: '
print 'fitted var: ',stdfit**2 ,'calculated var: ',variance


range = 5
x= np.linspace(-range,range,100)

n=0
h_poly = np.frompyfunc(mp.hermite,2,1)

const = 1/(np.sqrt(2**n * factorial(n)) * (2*np.pi)**0.25 )
rest = const* np.exp(-x**2/4)
tot  = (rest*h_poly(n,x/np.sqrt(2)))**2
plt.plot(x,tot,'r-')

#print np.trapz(tot,x)
sigma = 0.66
s = np.random.normal(0, sigma, 1000)
plt.hist(s,bins=50,normed=1)


#plt.plot(traj[:,0],traj[:,1])
#plt.plot(exact[:,0],exact[:,1],'r--')

#plt.plot(traj[:,0],np.exp(-gamma*traj[:,0]),'r--')

#x = mu + sigma*np.random.randn(1000)
#n, bins, patches = plt.hist(x, 50, facecolor='green', alpha=0.75)
#y = mlab.normpdf( bins, mu, sigma)
#l = plt.plot(bins, y, 'r--')


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



plt.show()
