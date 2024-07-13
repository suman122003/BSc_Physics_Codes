import numpy as np
from gaussn_prop import TDSE_time_evolution
from lho_fns import phi_n, pot_lho, psi_superposed

hbar, m = 1, 1
hbar_m = [hbar, m]

X = np.linspace(-5, 5, 100)
omega = 1
hbar_m = [hbar, m]
prpt = m, omega

def bndst_evolution(plot_):
    if plot_ == 'eigenstate':
        n = 5       # INPUT
        prphi = [hbar_m, omega, n]
        iters_solve, T_max, pause_time = 10, 6, 0.01
        TDSE_time_evolution(hbar_m, phi_n, prphi, pot_lho, prpt, X, 
                iters_solve, T_max, pause_time, method_='LU_solve_def', plot_=True)
    elif plot_ == 'superposed state':
        ns = [2,3,4,5,6,7]  # INPUT
        cs = [10, 6, 4, 5, 4, 7]   # INPUT
        prpsi = [hbar_m, omega, ns, cs]
        iters_solve, T_max, pause_time = 10, 20, 0.01
        TDSE_time_evolution(hbar_m, psi_superposed, prpsi, pot_lho, prpt, X, 
                iters_solve, T_max, pause_time, method_='LU_solve_def', plot_=True)

bndst_evolution(plot_='eigenstate')
# plot_: eigenstate, superposed state