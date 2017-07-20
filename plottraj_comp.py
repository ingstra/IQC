from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from numpy import pi as pi
from scipy.stats import norm
from scipy.misc import factorial


#traj_60 = np.loadtxt('traj_60.dat')
traj_no = np.loadtxt('traj.dat')
#traj_av40= np.loadtxt('traj_average40.dat')
#traj_av20= np.loadtxt('traj_average.dat')

traj = np.loadtxt('traj_mod_kappah.dat')

plt.plot(traj_no[:,0],traj_no[:,2],'k',linewidth=2,label=r'Modulate $\kappa_h$')
plt.plot(traj_no[:,0],traj_no[:,1],'k',linewidth=2)

plt.plot(traj[:,0],traj[:,2],'--r',linewidth=2)
plt.plot(traj[:,0],traj[:,1],'--r',linewidth=2,label=r'Modulate $n_h$')


#plt.plot(traj_no60[:,0],traj_no60[:,2],'--',linewidth=2,label=r'no noise 60')
#plt.plot(traj_no60[:,0],traj_no60[:,1],'--',linewidth=2,label=r'no noise 60')

#plt.plot(traj_no40[:,0],traj_no40[:,2],':',linewidth=2,label=r'no noise 40')
#plt.plot(traj_no40[:,0],traj_no40[:,1],':',linewidth=2,label=r'no noise 60')
#plt.plot(traj_60[:,0],traj_60[:,2],linewidth=2,label=r'60')
#plt.plot(traj_60[:,0],traj_60[:,1],linewidth=2,label=r'60')




#plt.plot(traj_av20[:,0],traj_av20[:,2],linewidth=2,label='average 20 trajs')
#plt.plot(traj_av40[:,0],traj_av40[:,2],linewidth=2,label='average 40 trajs')
#plt.plot(traj_60[:,0],traj_60[:,2],linewidth=2,label='average 60 trajs')

plt.legend()

plt.xlabel(r'ns')


plt.tight_layout()
#plt.savefig('testplot',figsize=(20,10))

plt.show()
