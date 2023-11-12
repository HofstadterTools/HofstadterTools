Peierls Phases
==============

In this section, we will demonstrate how the Peierls phases are computed for the conventional Hofstadter Model

.. math::
   H = -t\sum_{\braket{ij}} e^{i \theta_{ij}} c_i^\dagger c_j + \mathrm{H.c.},

where :math:`t` is the hopping amplitude, :math:`\braket{ij}` denotes nearest neighbors on a square lattice in the xy-plane, :math:`\theta_{ij}` is the Peierls phase acquired by hopping from site :math:`i` to :math:`j` due to a magnetic field in the z-direction, and :math:`c^{(\dagger)}` are the particle (creation)annihilation operators.

The Peierls phases acquired by hopping from site :math:`i\equiv(X_i,Y_i)` to site :math:`j\equiv(X_j,Y_j)` are given by

.. math::
   \theta_{ij} = \frac{2\pi}{\phi_0} \int_i^j \mathbf{A}\cdot \mathbf{dl},

where :math:`\phi_0=h/e` is the flux quantum, :math:`\mathbf{A}` is the vector potential, and :math:`\mathbf{dl}\equiv(\mathrm{d}x,\mathrm{d}y)` is an infinitesimal line element. There is gauge freedom in describing the magnetic field :math:`\mathbf{B}=\nabla\times\mathbf{A}`.

In HofstadterTools, in order to describe a perpendicular magnetic field in the z-direction :math:`\mathbf{B}=B\hat{\mathbf{e}}_z`, we choose the Landau gauge in the y-direction, i.e. with a conserved :math:`k_x` momentum, such that :math:`\mathbf{A}=-By\hat{\mathbf{e}}_x`. Parameterizing the line integral using

.. math::
   \begin{cases}
      x&=X_i + (X_j - X_i)\tau \\
      y&=Y_i + (Y_j - Y_i)\tau
   \end{cases}

where :math:`\tau\in[0,1]`, then yields

.. math::
   \int_i^j \mathbf{A}\cdot\mathbf{dl} = -\int_0^1 By (X_j-X_i) \mathrm{d}\tau.

Evaluating this integral yields

.. math::
   -\int_0^1 BY_i (X_j - X_i) \mathrm{d}\tau - \int_0^1 B (Y_j - Y_i)(X_j-X_i)\tau \mathrm{d}\tau = - BY_i(X_j-X_i) - \frac{B(Y_j - Y_i)(X_j - X_i)}{2},

hence, the Peierls phases may be written as

.. math::
   \theta_{ij} = -\left( \frac{2\pi B}{\phi_0} \right) \Delta X \left( Y_i + \frac{\Delta Y}{2} \right),

where :math:`\Delta Y = Y_j - Y_i` and :math:`\Delta X = X_j - X_i`. Note that, with this choice of gauge, the Peierls phases depend on *absolute* y coordinates but only *relative* x coordinates.

Since the magnetic flux density, :math:`n_{\phi}=BA/\phi_0`, is often defined with respect to the lattice unit cell area :math:`A=a^2`, and the lattice constant :math:`a` is often set to one, you may sometimes see this expression written as

.. math::
   \theta_{ij} = - 2\pi n_\phi \Delta X \left( Y_i + \frac{\Delta Y}{2} \right).

However, caution is required for this last step. Sometimes, especially when the area of a minimal hopping plaquette does not coincide with the lattice unit cell area, the flux density is defined with respect to a different area quantity.

This derivation of the Peierls phase may be readily generalized to the case of a general magnetic field :math:`\mathbf{B}=(B_x, B_y, B_z)` and a general gauge, e.g. Landau gauge in the x-direction :math:`\mathbf{A}=Bx\hat{\mathrm{e}}_y`, or symmetric gauge :math:`\mathbf{A}=(-By/2,Bx/2,0)`. :cite:`Peierls33`
