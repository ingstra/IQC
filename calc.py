import numpy as np
from scipy.integrate import quad
from scipy.misc import factorial
from scipy.special import eval_hermite as hermite
from scipy.stats import norm
import matplotlib.pylab as plt
import matplotlib.mlab as mlab


x=np.linspace(-5,5,100)

dx = x[1]-x[0]

def psi1(n,x):
	const = 1/(np.sqrt(2**n * factorial(n)) * (2*np.pi)**0.25 )
        rest = const* np.exp(-x**2/4)
        tot  = (rest*hermite(n,x/np.sqrt(2)))**2
        return tot

def psi2(n,x):
	return 1/(np.power(2,n)*factorial(n))*1/np.sqrt(np.pi)*np.exp(-np.power(x,2))*np.power(hermite(n,x),2);


n=0


psi1vec=psi1(n,x)
psi2vec=psi2(n,x)

print np.trapz(psi1vec)*dx,  np.trapz(psi2vec)*dx
print np.var(psi1vec,ddof=1),  np.var(psi2vec)

#plt.plot(x,psi1vec)
plt.plot(x,psi2vec)

pdf  = psi1vec

prob = pdf / pdf.sum() 
mu = x.dot(prob)
mom2 = np.power(x, 2).dot(prob)
var  = mom2 - mu**2
#print 'var 1', var

pdf  = psi2vec

prob = pdf / pdf.sum() 
mu = x.dot(prob)
mom2 = np.power(x, 2).dot(prob)
var  = mom2 - mu**2
print 'var 2', var

#################################3
current = np.loadtxt('current.dat')


values, bins, _ = plt.hist(current,50,normed=1,histtype='bar')

mufit, stdfit = norm.fit(current)
print 'fitted var: ',stdfit**2 


sigma = 0.66
s = np.random.normal(0, sigma, 1000)
#plt.hist(s,bins=50,normed=1)


plt.show()

#print 'psi2', quad(psi2,-np.inf,np.inf,args=n)



