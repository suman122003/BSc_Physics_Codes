import numpy as np
import numpy.linalg as npLA
import scipy.linalg as spLA
import matplotlib.pyplot as plt

def LU_solve_tridiag(A, b):
    '''
    Function to solve a matrix equation by LU factorization (efficient for tridiagonal matrices)
        ARGUMENTS
    A, b: matrices for for equation, AX = b
        RETURNS
    X: solution of the matrix equation
    '''
    A_arr, b_arr = np.array(A), np.array([b])
    ab = np.concatenate((A_arr, b_arr.T), axis=1) # augmented matrix
    n = len(ab)
    l = [[0 for j in range(n)] for i in range(n)]
    u = [[0 for j in range(n)] for i in range(n)]
    z = [0 for i in range(n)]
    for i in range(n-1):
        l[i][i-1] = ab[i][i-1]
        l[i][i] = ab[i][i] -l[i][i-1]*u[i-1][i]
        u[i][i+1] = ab[i][i+1]/l[i][i]
        z[i] = (ab[i][n]-l[i][i-1]*z[i-1])/l[i][i]
    l[n-1][n-2] = ab[n-1][n-2]
    l[n-1][n-1] = ab[n-1][n-1] -l[n-1][n-2]*u[n-2][n-1]
    z[n-1] = (ab[n-1][n]-l[n-1][n-2]*z[n-2])/l[n-1][n-1]
    X = [0 for i in range(n)]
    X[n-1] = z[n-1]
    for i in range(n-2, -1, -1):
        X[i] = z[i] - u[i][i+1]*X[i+1]
    return X

# SOLVING THE TDSE MATRIX EQUATION
def TDSE_mat_solve(hbar_m, psiR, psiI, V, dx, dt, method_):
    '''
        ARGUMENTS
    hbar_m: [hbar, m]
    psiR: (array) real part of wavefunction
    psiI: (array) imaginary part of wavefunction
    V: Potential
    dx, dt: steps in x, steps in t
    method_: numpy.linalg.inv, scipy.linalg.lu_solve
        RETURNS
    psiR: (array) real part of wavefunction
    psiI: (array) imaginary part of wavefunction
    '''
    hbar, m = hbar_m
    Nx = len(psiR)
    psicomp0 = [0 for i in range(Nx)]
    for i in range(Nx):
        psicomp0[i] = complex(psiR[i], psiI[i])
    al, bt = hbar*dt/(4*m*dx**2), dt/(2*hbar)
    gam = [2*al + bt*V[i] for i in range(Nx)]
    A = [[0 for j in range(Nx)] for i in range(Nx)]
    A[0][0], A[Nx-1][Nx-1] = 1, 1
    B = [[0 for j in range(Nx)] for i in range(Nx)]
    B[0][0], B[Nx-1][Nx-1] = 1, 1
    for j in range(1, Nx-1):
        A[j][j-1], A[j][j], A[j][j+1] = -al, gam[j]-1j, -al
        B[j][j-1], B[j][j], B[j][j+1] = al, -gam[j]-1j, al
    A, B, psicomp0 = np.array(A), np.array(B), np.array(psicomp0)
    b = B.dot(psicomp0)
    if method_ == 'numpy.linalg.inv':
        psicomp = (npLA.inv(A)).dot(b)
    elif method_ == 'numpy.linalg.solve':
        psicomp = npLA.solve(A, b)
    elif method_ == 'scipy.linalg.lu_solve':
        lu1, piv1 = spLA.lu_factor(A)
        psicomp = spLA.lu_solve((lu1, piv1), b)
    elif method_ == 'LU_solve_def':
        psicomp = LU_solve_tridiag(A, b)
    psiR = np.array([psicomp[i].real for i in range(Nx)])
    psiI = np.array([psicomp[i].imag for i in range(Nx)])
    return psiR, psiI

# NORMALIZATION
def psiNormRI(psiR, psiI, dx):    
    '''
        ARGUMENTS
    psiR: list or array of real psi values
    psiI: list or array of imaginary psi values
    dx: step in x
        RETURNS
    psiR: normalized list
    psiI: normalized list
    '''
    Nx = len(psiR)
    psiR, psiI = np.array(psiR), np.array(psiI)
    psimod2 = psiR**2 + psiI**2
    psiNorm1 = (np.sum(psimod2)*dx)**0.5
    psiR = [psiR[i]/psiNorm1 for i in range(Nx)]
    psiI = [psiI[i]/psiNorm1 for i in range(Nx)]
    return psiR, psiI

# ULTIMATE FUNCTION
def TDSE_time_evolution(hbar_m, wvfn, prwv, pot, prpt, X, iters_solve, T_max, pause_time, method_, return_=False, plot_=False):
    '''
        ARGUMENTS
    hbar_m: parameters [hbar, m]
    wvfn: initial wavefunction wvfn(prwv, X)
    prwv: parameters for wvfn
    pot: function for potential pot_fn(prpt, x)
    prpt: parameters for pot_fn
    X: the x array
    iters_solve: No. of iterations interval for real time plot
    T_max: max (final) time
    pause_time: Pause time plt.pause(tps)
    method_: numpy.linalg.inv, np.linalg.solve, scipy.linalg.lu_solve, LU_solve_def
    return_: True if the output is required
    plot_: True if plots are required
        RETURNS
    t_arr, E, psiR_arr, psiI_arr, psimd2_arr if return_=True
    None if return_ is not given
    (Wave function plot) if plot_=True
    '''
    hbar, m = hbar_m
    Nx = len(X)
    x0, xN = X[0], X[Nx-1]
    dx = (xN-x0)/(Nx-1)
    E, wvfnx = wvfn(prwv, X)
    V = pot(prpt, X)
    Vmx = max(max(V), abs(min(V)))
    dt = hbar/(hbar**2/(m*dx**2) + Vmx/2)   # not necessary
    psi = wvfnx     # initial
    psiR = [psi[i].real for i in range(Nx)]  # initial real part
    psiI = [psi[i].imag for i in range(Nx)]  # initial imaginary part
    psiR, psiI = psiNormRI(psiR, psiI, dx)
    t_arr, psiR_arr, psiI_arr, psimd2_arr = [], [], [], []
    Xmn, Xmx, Ymx = min(X), max(X), 1.5*max(np.abs(np.array(psiR)))

    # Rescaling
    if plot_ == True:
        if Vmx != 0:
            Efac = Ymx/(2*Vmx)
            Vplot = [V[i]*Efac for i in range(Nx)]
        print(f'Energy of the particle = {E}, Scaled Energy = {E*Efac}')
    t, itr, plti = 0, 1, 1
    while t <= T_max:
        if itr % iters_solve == 0:  # plotting at ts no. of iterations
            psimd2 = [psiR[i]**2+psiI[i]**2 for i in range(Nx)]
            if return_ == True:
                t_arr.append(t)
                psiR_arr.append(psiR), psiI_arr.append(psiI)
                psimd2_arr.append(psimd2)
            if plot_ == True:
                plt.axis([Xmn, Xmx, -Ymx, Ymx])
                if Efac != 0:
                    plt.plot(X, Vplot, ':k')
                    plt.axhline(E*Efac, label=f'E={E:.3}, E_plot={(E*Efac):.3}')
                    plt.fill_between(X, Vplot, facecolor='grey')
                plt.plot(X, psiR, label=r'$\psi_{real}(x,t)$')
                # print(f'Normalization check: {np.sum(psimd2)*dx}')
                plt.plot(X, psimd2, label=r'$|\psi(x,t)|^2$')
                plt.text(0.1*Xmx, -0.9*Ymx, f't = {t}')
                plt.legend(loc='upper right')
                plt.xlabel('$x$')
                # plt.savefig(f'plots/plot_{plti}')
                plt.pause(pause_time)
                plt.clf()
            plti += 1
        psiR, psiI = TDSE_mat_solve(hbar_m, psiR, psiI, V, dx, dt, method_)
        itr += 1
        t += dt
    if return_ == True:
        return np.array(t_arr), E, np.array(psiR_arr), np.array(psiI_arr), np.array(psimd2_arr)
    if plot_ == True:
        plt.show()
