import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')
theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
x = np.linspace(-20, 20, 1000)

y=x*0

z =  np.cos(x)
ax.plot(x, y, z, label='parametric curve')
#ax.legend()
plt.xlabel('x')
plt.ylabel('y')

#y2= np.sin(x)

#z2 = x*0
#ax.plot(x, y2, z2, label='parametric curve')
ax.fill(x)


plt.show()
