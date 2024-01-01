# Anharmonic Oscillator
# To find the eigen Value and eigen function of TISE By Numerov Method
# Shooting method
def Psi(mhdx2,psi,Vi,E):
    N=len(psi)
    psiE=[psi[i]   for i in range(N)]
    P=[mhdx2*(Vi[i]-E) for i in range(N)]

    for i in range(2,N):
        d=1-1/12*P[i]
        a=2*(1+5/12*P[i-1])
        b=-(1-1/12*P[i-2])
        psiE[i]=a/d*psiE[i-1]+b/d*psiE[i-2]
    return psiE


def numerovSchr(mhdx2,Vi,psi0,psi1,psiN,nodes,mxItr):
    N=len(Vi)-1
    Emx=max(Vi)
    Emn=min(Vi)
    psiIn=[0 for i in range(N+1)]
    psiIn[0],psiIn[1],psiIn[N]=psi0,psi1,psiN
    itr=0
    while abs(Emx-Emn) > 1e-6 and itr < mxItr:
        E=0.5*(Emx+Emn)
        psi=Psi(mhdx2,psiIn,Vi,E)

        # node Counting
        cnt=0
        for i in range(1,N-2):
            if psi[i]*psi[i+1] <0:
                cnt+=1
        if cnt > nodes:
            Emx=E
        elif cnt<nodes:
            Emn=E
        else:
            if psi[N-1]> psiN:
                Emn=E
            elif psi[N-1]< psiN:
                Emx=E
        itr+=1
    if itr <mxItr:
        return E,psi
    else:
        return None,None

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



def psiNorm(psi,dx):

    N=len(psi)
    psi2=[psi[i]**2 for i in range(N)]

    psiMod2=simp13xdis(dx,psi2)
    nrmPsi=[psi[i]/((psiMod2**0.5)*50000) for i in range(N)]
    return nrmPsi

#Solving eigen function of harmonic/anharmonic oscillator
# v(r)=1/2 k r^2 +1/3 b r^3

import matplotlib.pyplot as plt

def V(pr,r):
    k,b=pr
    return 0.5*k*r**2#+(1/3)*b*r**3

hbar,m=0.1,1
dr=0.005
mxItr=100
psi0,psiN=0,0
k=1
R0,RN=[0,0,0,0],[1.25,1.5,1.7,1.9]
psi0,psiN=0,0
for nodes in range(4):
    N=int((RN[nodes]-R0[nodes])/dr)
    dr=(RN[nodes]-R0[nodes])/N
    mhdx2=2*m*dr**2/hbar**2
    r=[R0[nodes]+i*dr for i in range(N+1)]
    Vi=[V([1,0.01],r[i]) for i in range(N+1)]
    psi1=((-1)**nodes)*1e-4
    E,psi=numerovSchr(mhdx2,Vi,psi0,psi1,psiN,nodes,mxItr)

    if E!=None:
        psi=psiNorm(psi,dr)

        plt.plot(r,psi,label=r'E=%.4f$\psi_%d(r)$'%(E,nodes))

xax=[0 for i in range(N+1)]
plt.plot(r,xax,'k')
plt.plot(r,Vi,'k--',label='Potential')
plt.legend(loc='best',prop={'size':12})
plt.xlabel('r')
plt.show()























    
