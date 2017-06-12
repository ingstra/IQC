from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from scipy.misc import factorial
import sympy.mpmath as mp
from mpl_toolkits.mplot3d import Axes3D

rho = np.loadtxt('rho.dat')

fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')

num_elements = len(rho)

xpos = [i for i in range(0,num_elements)]
ypos = [i for i in range(0,num_elements)]
xpos, ypos = np.meshgrid(xpos, ypos)
xpos = xpos.flatten()   # Convert positions to 1D array
ypos = ypos.flatten()
zpos = np.zeros(num_elements*num_elements)

dx = np.ones_like(zpos)
dy = np.ones_like(zpos)
dz = rho.flatten()

plt.xlabel('n')
plt.ylabel('m')

ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color='#00ceaa')


plt.show()
