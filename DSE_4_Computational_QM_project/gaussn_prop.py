import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as npLA
import scipy.linalg as spLA

from fn_TDSE_solve import TDSE_time_evolution

hbar, m = 1, 1

def gausswv(pr, X):
    x1, sig, k0 = pr
    E = (hbar**2/(2*m))*(k0**2 + 1/(2*sig**2))
    a = 1/((2*np.pi)**0.5 *sig)**0.5
    b = -1/(4*sig**2)
    Nx = len(X)
    psigs = [0 for xi in range(Nx)]
    for xi in range(Nx):
        gswv = a*np.exp(b*(X[xi]-x1)**2)
        if gswv < max(abs(X))*1e-10:
            gswv = 0
        gswv *= complex(np.cos(k0*X[xi]), np.sin(k0*X[xi]))
        psigs[xi] = gswv
    return E, psigs

def pot_free(V0, X):
    Nx = len(X)
    V = [V0 for i in range(Nx)]
    return V
def pot_step(pr, X):
    V0 = pr
    Nx = len(X)
    V = [0 for i in range(Nx)]
    V[int(Nx/2):] = [V0 for i in range(int(Nx/2))]
    return V
def pot_barrier(pr, X):
    V0, thk = pr
    Nx = len(X)
    V = [0 for i in range(Nx)]
    V[int(Nx/2):int(Nx/2)+thk] = [V0 for i in range(thk)]
    return V
def pot_periodic(pr, X):
    V0, b = pr
    Nx = len(X)
    a, n1 = 4*b, int(Nx/5)
    V = [0 for i in range(Nx)]
    while True:
        V[n1:n1+int(b)] = [V0 for i in range(int(b))]
        n1 += a+b
        if n1 + int(b) > Nx:
            break
    return V
def pot_nuclear(pr, X):
    Nx = len(X)
    Vd, Vc = pr
    V = [0 for i in range(Nx)]
    V[:int(Nx/3)] = [Vd for i in range(int(Nx/3))]
    V[int(Nx/3):] = [Vc*X[int(Nx/3)]/X[i] for i in range(int(Nx/3), Nx)]
    return V

def Gaussian_wavepacket_evolution(pot):
    if pot == 'pot_free':
        potfn = pot_free
        x0, xN, Nx = 0, 30, 100
        X = np.linspace(0, 30, 100)
        Nx = len(X)
        x0, xN = X[0], X[Nx-1]
        sig = 0.5
        x1, k0 = 5, np.pi
        V0 = 1
        prpt = V0
        pr = [hbar, m]
        prwv = [x1, sig, k0]
        iters_solve, T_max, pause_time = 2, 10, 0.01
        print('Free Particle plot')
    elif pot == 'pot_step':
        potfn = pot_step
        X = np.linspace(0, 50, 100)
        sig = 1.5
        x1, k0 = 15, 1.5*np.pi
        V0 = 5
        pr = [hbar, m]
        prwv = [x1, sig, k0]
        prpt = V0
        iters_solve, T_max, pause_time = 3, 30, 0.01
        print('Potential Step plot')
    elif pot == 'pot_barrier':
        potfn = pot_barrier
        X = np.linspace(0, 40, 100)
        sig = 1
        x1, k0 = 12, np.pi
        V0, thk = 5, 10
        pr = [hbar, m]
        prwv = [x1, sig, k0]
        prpt = [V0, thk]
        iters_solve, T_max, pause_time = 3, 21, 0.01
        print('Potential Barrier plot')
    TDSE_time_evolution(pr, gausswv, prwv, potfn, prpt, X, 
                iters_solve, T_max, pause_time, method_='LU_solve_def', plot_=True)
    return None
