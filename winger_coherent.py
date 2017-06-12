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


x=np.linspace(-4,4,500)
p=x
rho=fock(4,0)
W=wigner(rho,x,p)
X, P = np.meshgrid(x,p)

def c(x,p):
    return np.exp(-((x-2)**2 + (p-2)**2))/np.pi
    
wcmap=w_cmap(W,shift=0)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, P, c(X,P),cmap=wcmap)
plt.tight_layout()

fig.colorbar(surf, shrink=0.7)

plt.xlabel(r'$x$')
plt.ylabel(r'$p$')

range= [-6,6]
print integrate.nquad(c,[range,range])

plt.savefig('figures/coherent.pdf',figsize=(10,10))

plt.show()



