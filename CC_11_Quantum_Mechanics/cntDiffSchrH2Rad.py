##Radial Wave Function of Hydrogen Atom
# D2y/dr2= -2/r* dy/dr-2m/h2*(E-v(r))y
# v(r)=e2/4 pi E0  (E0=epsilon0)




# Simpson's 1/3 rd rule for discrete function
def simp13xdis(h,fx):
    n=len(fx)
    I=0
    for i in range(n):
        if i==0 or i==n-1: # i 1st or last
            I+=fx[i]
        elif i%2 !=0: # i odd
            I+=4*fx[i]
        else: #i even
            I+=2*fx[i]
        I=I*h/3
    return I



#  Normalization of discrete wavefunction
#from  import simp13Xdis
def psiNorm(psi,dx):
    import numpy as np
    N=len(psi)
    psi2=[psi[i]**2 for i in range(N)]
    psiMod2=simp13xdis(dx,psi2)
    nrmPsi=[psi[i]/(np.sqrt(psiMod2)*160000) for i in range(N)]
    return nrmPsi



# Propagator with Central Difference Method Calculation
def propCntDiff(pr,p,q,r,x,y,dx):
    N=len(x)
    yy=[y[i] for i in range(N)]
    for i in range(2,N):
       a=2+dx**2*q(pr,x[i-1])
       b=-(1+dx/2.0*p(x[i-1]))
       c=(dx**2)*r(x[i-1])
       d=1-dx/2*p(x[i-1])
       yy[i]=a/d*yy[i-1]+b/d*yy[i-2]+c/d
    return yy


def CntDiffEigVal(prMn,prMx,p,q,r,x0,y0,xN,yN,y1,dx,nodes,mxItr):
        N=int((xN-x0)/dx)
        dx=(xN-x0)/N
        x=[x0+i*dx for i in range(N+1)]
        y=[0 for i in range(N+1)]
        y[0],y[1],y[N]=y0,y1,yN
        itr=0
        while abs(prMx-prMn) > 1e-6 and itr < mxItr:
            pr=0.5*(prMx+prMn)
            yy=propCntDiff(pr,p,q,r,x,y,dx)
            cnt=0
            for i in range(1,N-2):
                if yy[i]*yy[i+1] <0:
                    cnt+=1
            if cnt > nodes:
                prMx=pr
            elif cnt < nodes:
                prMn=pr
            else:
                if yy[N-1] >yN:
                 prMn=pr
                elif yy[N-1] < yN:
                 prMx=pr

                 itr+=1
        if itr < mxItr :
            return pr,x,yy
        else:
            return None, None,None

from math import *
import matplotlib.pyplot as plt

hbar,m=0.1,1
e2=0.2

mh2=2*m/hbar**2
# potential
def V(e2,r):
    return  (exp(-6*r)-exp(-3*r)) # -e2/r
def P(r):
    return -2/r
def Q(E,r):
    return -mh2*(E-V(e2,r))
def R(r):
    return 0
stln=['k','k--','k:','k-.']
dr=0.01
mxItr=100
r0,psi0,rN,psiN=1e-6,0,[0.7,1.2,2.2,3.2],0
# First four eigenfunctions
for nodes in range(4):
    Emn,Emx=V(e2,r0),V(e2,rN[nodes])
    psi1=psi0+((-1)**nodes)*1e-4
    E,r,psi=CntDiffEigVal(Emn,Emx,P,Q,R,r0,psi0,rN[nodes],
                          psiN,psi1,dr,nodes,mxItr)

            
    if E!=None:  # if solution is obtained
        psi=psiNorm(psi,dr)

        rPsi=[r[i]*psi[i] for i in range(len(r))]
        plt.plot(r,rPsi,stln[nodes],label=r'$r\psi_%d(r)$ for $E_%d$ = %.3f' %
                 (nodes+1,nodes+1,E))
        plt.xlabel('r')
        plt.legend(loc='best',prop={'size':12})
plt.show()




        
