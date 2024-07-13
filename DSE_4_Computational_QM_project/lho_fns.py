import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

def phi_n(pr, X):
    '''
        ARGUMENTS
    pr: [hbar, m]
    X: X array
        RETURNS
    En: Energy of eigenstate
    eigf: eigenstate as an array
    '''
    hbar_m, omega, n = pr
    hbar, m = hbar_m
    Vx = 0.5*m*(omega**2)* X**2
    Nx = len(X)
    dx = X[1] - X[0]
    D_mat = (np.diag(-2*np.ones(Nx)) + np.diag(np.ones(Nx-1), 1)
                        + np.diag(np.ones(Nx-1), -1))/dx**2
    T_mat = (-hbar**2/(2*m)) * D_mat
    V_mat = np.diag(Vx*np.ones(Nx))
    H_mat = T_mat + V_mat
    eigenvals1, eigenvecsT1 = eigh(H_mat)
    En = eigenvals1[n]
    eigf = eigenvecsT1[:, n]
    eigf = eigf/np.sum(np.abs(eigf)**2*dx)**0.5
    return En, eigf

def pot_lho(pr, X):
    m, omega = pr
    return (1/2)*m*omega*np.array(X)**2

def psi_superposed(pr, X):
    '''
        ARGUMENTS
    pr: [[hbar, m], omega, n_vals, c_vals] for psi = sum(c*phi_n)
    X: X array
        RETURNS
    Es: Energy of superposed state
    psis: Superposed state as an array
    '''
    hbar_m, omega, n_vals, c_vals = pr
    ns, cs = n_vals, np.array(c_vals)
    cs = cs/np.sum(cs**2)**0.5
    Es = 0
    psis = np.array([0 for i in range(len(X))], dtype=float)
    for i in range(len(ns)):
        prphi = hbar_m, omega, ns[i]
        Ei, psii = phi_n(prphi, X)
        Es += (cs[i])**2 *Ei
        psii = cs[i]*psii
        psis += psii
    return Es, psis
