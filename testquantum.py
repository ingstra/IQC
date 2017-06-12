import numpy as np
from qutip import *
import matplotlib as mlp
from matplotlib import cm
import matplotlib.pylab as plt

rho= fock_dm(2)
print(rho)

xvec = np.linspace(-5,5,200)

W = wigner(rho, xvec, xvec)


rhom = [[1,0],[0,0]]
print rhom

rhomq = Qobj(rhom)
print rhomq
