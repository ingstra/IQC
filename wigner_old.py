from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from matplotlib import cm
from scipy.stats import norm
from scipy.misc import factorial
from scipy.special import eval_genlaguerre as L
from scipy.special import eval_laguerre as lag
from scipy import integrate
import sympy.mpmath as mp
from mpl_toolkits.mplot3d import Axes3D

rho = np.loadtxt('rho.dat')

rho=np.array([[0, 0],[0,1]])
NFock = len(rho)

x=np.linspace(-6,6,500)
p=x

x, p = np.meshgrid(x,p)

laguerre = np.vectorize(L)

def W_mn(x,p,m,n):
    if m >= n:
        w = np.exp(-x*x-p*p)*np.power(-x+1j*p,m-n)*np.sqrt(np.power(2,m-n)*factorial(n)/factorial(m))* np.power(-1,m)*laguerre(n,m-n,2*x*x+2*p*p)/np.pi
    elif m < n: 
        w = np.exp(-x*x-p*p)*np.power(x+1j*p,n-m)*np.sqrt(np.power(2,n-m)*factorial(m)/factorial(n))* np.power(-1,n)*laguerre(n,m-n,2*x*x+2*p*p)/np.pi
    return w   

def W_mn2(x,p,m,n):
     if m >= n:
        w= np.power(-1,n)*np.sqrt(factorial(n)/factorial(m)) *np.exp(-x*x-p*p)*np.power(2*p*p+2*x*x,(m-n)/2)*laguerre(n,m-n,2*x*x+2*p*p)*np.exp(1j*(m-n)*np.arctan(p/x)) /np.pi
     elif m < n: 
        w= np.power(-1,m)*np.sqrt(factorial(m)/factorial(n)) *np.exp(-x*x-p*p)*np.power(2*p*p+2*x*x,(n-m)/2)*laguerre(n,n-m,2*x*x+2*p*p)*np.exp(1j*(n-m)*np.arctan(p/x))/np.pi
     return w 

def wigner(x,p):
    W=0
    for m in xrange(0,NFock):
        for n in xrange(0,NFock):
            W = W + rho[m][n]*W_mn(x,p,m,n)     
            
    return np.real(W)



def wigner2(x,p):
    W=0
    for m in xrange(0,NFock):
        for n in xrange(0,NFock):
            W = W + rho[m][n]*W_mn2(x,p,m,n)     
    return W       

lag2=np.vectorize(lag)

def W_n(x,p,n):
    return 2*np.power(-1,n)*np.exp(-x*x-p*p)*lag(n,2*x*x+2*p*p)/np.pi

print rho

W1=wigner(x,p)

W2=wigner2(x,p)



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(x, p, np.real(W2),cmap=cm.coolwarm)
#ax.plot_surface(x, p, np.real(W2),cmap=cm.coolwarm)


#print integrate.nquad(wigner2,[[-6,6],[-6,6]])
#print integrate.dblquad(wigner,-6,6, lambda x: -6, lambda x: 6)


plt.show()


