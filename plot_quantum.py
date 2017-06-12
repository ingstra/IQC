import numpy as np
from qutip import *
import matplotlib as mlp
from matplotlib import cm
import matplotlib.pylab as plt

rho= fock_dm(2)
print(rho)

xvec = np.linspace(-5,5,200)

W = wigner(rho, xvec, xvec)

#fig = plt.figure()

#ax1 = fig.add_subplot(111)

#cont0 = ax.contourf(xvec, xvec, W, 100)
#----------------------------------------
#wmap = wigner_cmap(W) 

#nrm = mlp.colors.Normalize(-W.max(), W.max())

#fig, axes = plt.subplots(1, 2, figsize=(10, 4))

#plt1 = axes[0].contourf(xvec, xvec, W, 100, cmap=cm.RdBu, norm=nrm)

#axes[0].set_title("Standard Colormap");

#cb1 = fig.colorbar(plt1, ax=axes[0])

#plt2 = axes[1].contourf(xvec, xvec, W, 100, cmap=wmap)  # Apply Wigner colormap

#axes[1].set_title("Wigner Colormap");

#cb2 = fig.colorbar(plt2, ax=axes[1])

#------------------------
N = 2

lbls_list = [[str(d) for d in range(N)], [r'$\ket{a}$', "d"]]
xlabels = []
for inds in tomography._index_permutations([len(lbls) for lbls in lbls_list]):
     xlabels.append("".join([lbls_list[k][inds[k]]
                for k in range(len(lbls_list))]))
     
                                   
xlabels = [r'$\bra{0}$',r'$\bra{1}$'] 
ylabels = [r'$\ket{0}$',r'$\ket{1}$'] 
fig, ax = matrix_histogram(rho, xlabels, ylabels, limits=[-2,2])
ax.view_init(azim=-55, elev=45)

#rho_ss = steadystate(rho, [np.sqrt(0.1) * a, np.sqrt(0.4) * b.dag()])



fig, ax = hinton(rho) # xlabels=xlabels, ylabels=xlabels)






plt.show()
