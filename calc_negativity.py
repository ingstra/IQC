from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import matplotlib as mlp
import numpy as np
from matplotlib import cm
from scipy.stats import norm
from scipy.misc import factorial
from scipy.special import eval_genlaguerre as L
from scipy import integrate

T=[0.2*i for i in range (1,51)]
Omega=[0.1*i for i in range (1,51)]

x=np.linspace(-3,3,500)
p=x
        
X, P = np.meshgrid(x,p)
dx=x[1]-x[0]
dp=dx
laguerre = np.vectorize(L)
        
def W_mn(x,p,m,n):
     if n>=m:
         w = np.exp(-x*x-p*p)* np.power(-1,m)*np.sqrt(np.power(2,n-m)*factorial(m)/factorial(n))*np.power(x-1j*p,n-m)*laguerre(m,n-m,2*x*x+2*p*p)/np.pi
     elif n < m: 
         w = np.exp(-x*x-p*p)* np.power(-1,n)*np.sqrt(np.power(2,m-n)*factorial(n)/factorial(m))*np.power(x+1j*p,m-n)*laguerre(n,m-n,2*x*x+2*p*p)/np.pi
     return w   

N = [[0 for x in range(1,51)] for y in range(1,51)] 

print len(N)
for n in T:
    for m in Omega:
     #   filename_im = "run1/rho_im_T" + str(n) + "_Omega" + str(m)
      #  filename_re = "run1/rho_re_T" + str(n) + "_Omega" + str(m)

       # rho_real = np.loadtxt(filename_re)
        #rho_im = np.loadtxt(filename_im)

#        rho = rho_real + 1j*rho_im

 #       NFock = len(rho)
          
        def wigner_calc(x,p):
            W=0
            for m in xrange(0,NFock):
                for n in xrange(0,NFock):
                    W = W + rho[m][n]*W_mn(x,p,m,n)     
                    
                    return np.real(W)

                
        def wigner_calc_abs(x,p):
            W=0
            for m in xrange(0,NFock):
                for n in xrange(0,NFock):
                    W = W + rho[m][n]*W_mn(x,p,m,n)     
                    
                    return np.abs(np.real(W))
                
        def wignerint(x,p):
            return wigner_calc_abs(x,p) - wigner_calc(x,p)
                                                       
                            
        W=wigner_calc(X,P)


        Wabs=np.abs(W)
        int1 = np.trapz(Wabs-W)*dp
        int2=np.trapz(int1)*dx
        print int2
        N[int(n/0.2)-1][int(m/0.1)-1] = int2

print N
