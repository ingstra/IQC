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
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from qutip import *

from wigner_cmap import w_cmap



rho_real = np.loadtxt('rho_real.dat')
rho_im = np.loadtxt('rho_im.dat')

rho = rho_real + 1j*rho_im

rho=np.array([[0, 0],[0,1]])
NFock = len(rho)

x=np.linspace(-3,3,500)
p=x

X, P = np.meshgrid(x,p)

laguerre = np.vectorize(L)

def W_mn(x,p,m,n):
    if n>=m:
          w = np.exp(-x*x-p*p)* np.power(-1,m)*np.sqrt(np.power(2,n-m)*factorial(m)/factorial(n))*np.power(x-1j*p,n-m)*laguerre(m,n-m,2*x*x+2*p*p)/np.pi
    elif n < m: 
        w = np.exp(-x*x-p*p)* np.power(-1,n)*np.sqrt(np.power(2,m-n)*factorial(n)/factorial(m))*np.power(x+1j*p,m-n)*laguerre(n,m-n,2*x*x+2*p*p)/np.pi
    return w   

def wigner_calc(x,p):
    W=0
    for m in xrange(0,NFock):
        for n in xrange(0,NFock):
            W = W + rho[m][n]*W_mn(x,p,m,n)     
            
    return np.real(W)


print rho

W=wigner_calc(X,P)

wcmap=w_cmap(W,shift=-1e-2)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, P, np.real(W),cmap=wcmap)
plt.tight_layout()

#fig.colorbar(surf, shrink=0.7)

plt.xlabel(r'$x$')
plt.ylabel(r'$p$')



plt.show()


