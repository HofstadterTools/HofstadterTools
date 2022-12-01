Peierls Phases
==============

Using a simple example, we will demonstrate how the Peierls phases are computed for the Hofstadter Model

.. math::
   H = \sum_{\braket{ij}} e^{i \theta_{ij}} c_i^\dagger c_j + \mathrm{H.c.},

where :math:`\braket{ij}` denotes nearest neighbors on a square lattice in the xy-plane, :math:`t` are the hopping amplitudes, :math:`\theta_{ij}` are the Peierls phases, and :math:`c^{(\dagger)}` are the (creation)annihilation operators for spinless fermions :cite:`Andrews20`.

The Peierls phases acquired by hopping from site :math:`i\equiv(X_i,Y_i)` to site :math:`j\equiv(X_j,Y_j)` are given by

.. math::
   \theta_{ij} = \frac{2\pi}{\phi_0} \int_i^j \mathbf{A}\cdot \mathbf{dl},

where :math:`\phi_0=h/e` is the flux quantum, :math:`\mathbf{A}` is the vector potential, and :math:`\mathbf{dl}\equiv(\mathrm{d}x,\mathrm{d}y)` is an infinitesimal line element. There is gauge freedom in describing the magnetic field :math:`\mathbf{B}=\nabla\times\mathbf{A}`.

For example, in order to describe a perpendicular magnetic field in the z-direction :math:`\mathbf{B}=B\hat{\mathbf{e}}_z`, we may choose the Landau gauge in the x-direction :math:`\mathbf{A}=Bx\hat{\mathbf{e}}_y`. Parameterizing the line integral using

.. math::
   \begin{cases}
      x&=X_i + (X_j - X_i)\tau \\
      y&=Y_i + (Y_j - Y_i)\tau
   \end{cases}

where :math:`\tau\in[0,1]`, yields

.. math::
   \int_i^j \mathbf{A}\cdot\mathbf{dl} = \int_0^1 Bx (Y_j-Y_i) \mathrm{d}\tau.

Evaluating this integral yields

.. math::
   \int_0^1 BX_i (Y_j - Y_i) \mathrm{d}\tau + \int_0^1 B (X_j - X_i)(Y_j-Y_i)\tau \mathrm{d}\tau = BX_i(Y_j-Y_i) + \frac{B(X_j - X_i)(Y_j - Y_i)}{2},

hence, the Peierls phases may be written as

.. math::
   \theta_{ij} = \left( \frac{2\pi B}{\phi_0} \right) (Y_j - Y_i) \left( X_i + \frac{X_j - X_i}{2} \right).

Note that, with this choice of gauge, the Peierls phases depend on *absolute* x coordinates but only *relative* y coordinates.

Moreover, since the magnetic flux per unit cell (flux density) is given as :math:`n_{\phi}=BA/\phi_0`, where :math:`A=a^2` is the unit cell area with lattice constant :math:`a=1`, we may write

.. math::
   \theta_{ij} = 2\pi n_\phi (Y_j - Y_i) \left( X_i + \frac{X_j - X_i}{2} \right).

This expression may be readily generalized to the case of a general magnetic field :math:`\mathbf{B}=(B_x, B_y, B_z)` and a general gauge, e.g. Landau gauge in the y-direction :math:`\mathbf{A}=-By\hat{\mathrm{e}}_x`, or symmetric gauge :math:`\mathbf{A}=(-By/2,Bx/2,0)`.
