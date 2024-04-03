# CC 05 - Mathematical Physics 2


## Sem 3 all codes in syllabus [file](sem3_all_codes.ipynb)
(ipynb and html file)

> **Bessel Functions** (by available function, calculations and recurrence formula)

> **Legendre Polynomials** (by available function, calculations and recurrence formula)

> **Least Square Method** (for fitting a straight line)

> **Gauss Elimination Method**

> **Inverse of a matrix by Gauss elimination method**

> **Solving Differential Equations by Euler's method**

> **Modified Euler Method for solving ODE**

> **Euler's metod for 2nd order ODE** (example - damped harmonic oscillation)

## Polynomial fit

Here a function is defined to fit a polynomial in some given data.

For a fitting polynomial, 
$$y = a_0 +a_1 x +a_2 x^2 +... +a_n x^n = \sum_{i=0}^n a_i x^i$$
The fitting parameters can be calculated by the matrix relation (for $n=3$);
$$\begin{bmatrix}
a_0 \\
a_1 \\
a_2 \\
\end{bmatrix}
=
\begin{bmatrix}
n & \sum x & \sum x^2 \\
\sum x & \sum x^2 & \sum x^3 \\
\sum x^2 & \sum x^3 & \sum x^4 \\
\end{bmatrix} ^{-1}
\begin{bmatrix}
\sum y \\
\sum xy \\
\sum x^2y \\
\end{bmatrix}$$

This idea can be extended for higher degree polynomials. We can do this by available functions, `polyfit` (in `numpy`) and `curve_fit` (in `scipy.optimize`) also.


## Special Functions [file](Special_Functions_SKP.ipynb)

Some special functions available in SciPy and SymPy are given here.

* **Gamma function**
* **Beta function**
* **Error function and complementary error function**

* **Legendre Polynomials**
* **Bessel functions**
* **Hermite Polynomials**
* **Laguerre Polynomials**

* **Permutations and Combinations**
* **Riemann Zeta function**

