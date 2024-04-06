# CC 11 - Quantum Mechanics

## 5th Semester all codes in Syllabus
(.ipynb and .html files)

> **Eigenvlaue Problems**

> **Central Difference Method** (Central potential, Potential well: infinite and finite, triangular potential, Radial equation: H-atom)

> **Numerov Method** (Infinite potential well, Linear Harmonic Oscillator, LHO with perturbation)

> **Exam Questions - 2023** (The questions are solved using Central difference method. Solution of all the 5 questions are given here.)



## DB sir - codes based on lectures [file](QM_DB_SKP.ipynb)

All the functions are created by me. Here you can explore codes of the following topics:

* L-S Coupling: Get the results of LS coupling.
* J-J Coupling: Get the results of JJ coupling.

* Selection rules: Here the selection rules are applied for a transition.
* Anamolous Zeemann Effect: Get energy levels.

## QM IITD - Dirac Delta Function [file](QM_iitd_1_Dirac_delta_function_skp.ipynb)

## QM IITD - Harmonic Oscillator and Spherical Harmonics [file](QM_iitd_2.ipynb)

Codes for following topics are given here:
* Harmonic Oscillator
* Spherical Harmonics

## QM IITD - Theory of Angular Momentum [file](QM_iitd_3_angular_momentum_skp.ipynb)

Code for following topics are given here:
* Pauli spin matrices
* Angular momentum operators for a given j
* Addition of angular momentum - Clebsch Gordan coefficients

## Time evolution of Gaussian wavepacket [file](time_evolution_skp.ipynb)

Time evolution of Gaussian wavepacket is shown here. You can read the theory on the book *Quantum Mechanics by Zettili* (section 1.8). Here all the steps are done for the time evotion of wavefunction for a free particle. The Fourier transformations and inverse Fourier transformations are done by integrations in this code. At the end, you can notice the broadening and deformation of wavefunction (as well as probability density) with the increase of time.

## H-atom problem - Symbolic Solution and visualization [file](H_atom_sp_visualization_skp.ipynb)

Here the solutions of the Hydrogen atom problem is given. Explore the solutions and graphs.

- I have defined the functions in the **Defining Function** [section](#defining-function). Here I have taken functions from `sympy.physics.hydrogen` and created numerical functions from them. The defined function, named **`hydrogen_wavefn_all`** gives the analytical and graphical view of total wavefunction $\psi_{nlm}(r, \theta, \phi)$ and radial wavefunction $R_{nl}(r)$ and **`hydrogen_wavefn_return`** returns 2 functions  **`phi_nlm_fn(r, theta, phi)`** and   **`R_nl_fn(r)`** for those wavefunctions. 

- In the **Give Inputs and Get States** [section](#give-inputs-and-get-states), you can give arbitary values of (n, l, m) and get the states. Here you need to set the number of grid points (`N1`) and the maximum limit of radius `r1` to view the total graphs. 

- Once you are done with the above things, go to the **Working with Functions** [section](#working-with-functions). Here I have used the returned functions (**`phi_nlm_fn(r, theta, phi)`** and   **`R_nl_fn(r)`**) and done few more plots from these.

- At last, in the **Effects of n, l and m** [section](#effects-of-n-l-and-m), the variations of n, l and m are shown. Here you can see how wavefunctions change (analytically and graphically) when n, l and m changes.

## Angular momentum matrices [file](ang_mom_matrices_skp.ipynb)

The codes are done to calculate angular momentum matrices using `sympy`.

* Defining functions
  - $J^2 \, \implies \,$ **`J2_mat(j)`**
  - $J_z \, \implies \,$ **`Jz_mat(j)`**
  - $J_x \, \implies \,$ **`Jx_mat(j)`**
  - $J_y \, \implies \,$ **`Jy_mat(j)`**
  - $J_+ \, \implies \,$ **`J_plus_mat(j)`**
  - $J_- \, \implies \,$ **`J_minus_mat(j)`**
  
* Angular momentum - matrix form (input total angular momentum and get the matrices)
* Eigenvalues and Eigenvectors
(Output - $ \left[ \,\left( \ eigenvalue, \  degeneracy, \  \left[ \, eigenvectors \ \right] \,\right) \,\right] $)

## Quantum Mechanics - 1D analysis [file](QM1d_analysis_SKP.ipynb)

* At first, the [Theory](#theory) behind this code is given here.

* Different potentials are defined in the Potentials [section](#potentials). Choose any 1 of the potentials.
* Then skip to the Plotting [section](#plotting-wavefunction-and-probability-density). Here, by watching the plots for near about 15-20 eigenstates, you'll have some interpretetion about the effect of the potential. Now, select some of the states of which you want to see the superposition. To do so, go to the Superposition [section](#superposition). 
* Now, you can also get the wavefunction in momentum space in the Momentum space [section](#momentum-space). It's a tedious analytical task to calculate expectation values and uncertainty for a given potential. You can do that in seconds in the Expectation values and Uncertainties [section](#expectation-values-and-uncertainty) (you can use the potential for 1d linear harmonic oscillator and get the uncertainty and it's totally correct, indicating the correctness of this code). 
* If you are only interested in the plots and calculations then that's all in this notebook. Remember, wherever you are needed to give some input, I have used ***`# INPUT`*** in that line.
* The Time evolution [section](#time-evolution) is not developed, so skip it. There is a module named `qutip` for the works related to Quantum Mechanics. I have used that for finding eigenstates from the hamiltonian in the QuTip analysis [section](#qutip-analysis). But it's suggested to skip this section. 
* At last, in the Roughs [section](#roughs), I have discussed the comparison of different python functions by showing their accuracy and runtime.



