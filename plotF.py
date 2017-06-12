from __future__ import division

import matplotlib.pylab as plt
import numpy as np
from scipy.fftpack import fft
from scipy.integrate import quad



z=np.linspace(-10,10,10000)
w=z
nr=1000000
x= np.linspace(-25,25,nr)

T=3



Omega = 0
gamma=1
gammafilt=6
w0=0

def F(w):
    return (1j-1j*np.exp(T*1j*w))/(np.sqrt(2*np.pi)*w)

filter = F(z)


def step(x):
   return 1 * (0 < x)*(x < 0.5)

def Lorentz(w):
    return gammafilt**2/(gammafilt**2+(w)**2)

def Lorentz2(w):
    return gamma**2/(gamma**2+(w)**2)

L=Lorentz(w)
L2 = Lorentz2(w)

def S(w):
    term1 = (1/4)*gamma/((gamma/2)**2+w**2)
    term2 = (3/16)*gamma/((3*gamma/4)**2+(w-Omega)**2)
    term3 = (3/16)*gamma/((3*gamma/4)**2+(w+Omega)**2)
    return term1+term2+term3

def minim(x):
    return min(S(x),F(x))

def const(x):
    return 1

spectrum=S(w)
lmin=np.minimum(spectrum,L)
l2min=np.minimum(spectrum,L2)
fmin=np.minimum(spectrum,filter)
print 'lorentz overlap: ',np.trapz(lmin),'sinc overlap: ',np.trapz(fmin), 'wrong filter overlap', np.trapz(l2min)

    

blah = step(x)
bleh = fft(blah)

plt.plot(w,S(w),label='spectrum')
plt.plot(w,L,label='wrong $\gamma$')
plt.plot(w,L2,label='$\gamma=1$')

#plt.plot(w,filter,label='filter')
plt.xlabel(r'$\omega$')
plt.ylabel(r'$S(\omega)$')

plt.legend()
plt.tight_layout()


plt.savefig('./figures/filter2_gamma.pdf',figsize=(10,10))
plt.show()
