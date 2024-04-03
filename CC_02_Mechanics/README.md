# CC 02 - Mechanics

## Moment of Inertia Tensor (discrete system) [file](MOI_tensor_skp.ipynb)

The angular momentum is
$$ L = I \omega $$
The Moment of Inertia tensor is;
$$ I =\begin{bmatrix}
I_{xx} & I_{xy} & I_{xz} \\
I_{yx} & I_{yy} & I_{yz} \\
I_{zx} & I_{zy} & I_{zz} \\
\end{bmatrix}$$
where,
$$ I_{xx} = \sum_i m_i (r_i^2 - x_i^2)$$
$$ I_{yy} = \sum_i m_i (r_i^2 - y_i^2)$$
$$ I_{zz} = \sum_i m_i (r_i^2 - z_i^2)$$
$$ I_{xy} = I_{yx} = -\sum_i m_ix_iy_i $$
$$ I_{yz} = I_{zy} = -\sum_i m_iy_iz_i $$
$$ I_{zx} = I_{xz} = -\sum_i m_iz_ix_i $$
$$(r_i = \sqrt{x_i^2 + y_i^2 + z_i^2})$$

By using this concept, we can calculate MOI of discrete mass system very easily. To calculate this for continuous mass (rod, disc etc.), the challenging task is to create a discrete distribution that can be able to represent a continuous mass distribution.

Here, moment of inretia is calculated for
* Discrete masses
* Linear mass (rod)
* Annular disc
* Disc

