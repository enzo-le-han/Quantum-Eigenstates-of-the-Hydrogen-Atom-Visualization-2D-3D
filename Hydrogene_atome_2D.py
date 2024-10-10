import numpy as np
from scipy.special import factorial, lpmn
import matplotlib.pyplot as plt


def Laguerre(alpha,degree,x): 
    if degree == 0:
        return 1 
    if degree == 1:
        return -x+alpha+1   
    else:
        L = 0
        L0, L1 = 1, -x+alpha+1
        
        for i in range(2,degree+1):   
            L = (2+(alpha-1-x)/i)*L1 - (1+(alpha-1)/i)*L0
            L0, L1 = L1, L
        return L

def Legendre(m,l,x):
    return lpmn(m,l,x)[0][m][l]   
Legendre_ = np.vectorize(Legendre, excluded=['m', 'l'])


def Radial(r,n,l):
    rho = 2*r/n
    return (2/n**2)*np.sqrt(factorial(n-l-1)/factorial(n+l))*(rho**l)*Laguerre(2*l+1,n-l-1,rho)*np.exp(-rho/2)

def Ylm(t,l,m):
    if m >= 0:
        N = np.sqrt((2*l+1)*factorial((l-m))/(4*np.pi*factorial((l+m))))
        return N * Legendre_(m,l,np.cos(t)) * (-1)**m
    else:
        N = np.sqrt((2*l+1)*factorial((l+m))/(4*np.pi*factorial((l-m))))
        return N * Legendre_(-m,l,np.cos(t))    



length = 40
pixels = 200
x = np.linspace(-length,length,pixels)
z = np.linspace(-length,length,pixels)

X, Z = np.meshgrid(x, z)
R = np.sqrt(X**2+Z**2)
A = np.arccos(Z/R)

n = 4
l = 2
for m in range(-l,l+1):
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    F = (Radial(R,n,l)*Ylm(A,l,m))**2   

    max = np.max(F)
    pc = ax.imshow(F,cmap='magma',vmax=max/3)
    ax.set_title("n,l,m = "+str(n)+', '+str(l)+', '+str(m))
    ax.axis('off')
    fig.tight_layout()
    plt.show()


# Radial function
R_fct = []  
Proba = [] 
r = np.linspace(0,35,1000)

for n in (1,2,3):
    for l in range(0,n):
        fig, axs = plt.subplots(1, 2,figsize=(10, 4))
        R_fct.append(Radial(r,n,l))
        Proba.append(Radial(r,n,l)**2*r**2)

        axs[0].plot(r, R_fct[-1])
        axs[0].set_title("Radial function n,l = "+str(n)+', '+str(l))
        axs[1].plot(r, Proba[-1])
        axs[1].set_title("Radial probability density n,l = "+str(n)+', '+str(l))

        fig.tight_layout()
        plt.show()