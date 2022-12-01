Band Structures
===============

Using a simple example, we will demonstrate how the band structures are computed for the Hofstadter Model

.. math::
   H = \sum_{\braket{ij}} e^{i \theta_{ij}} c_i^\dagger c_j + \mathrm{H.c.},

where :math:`\braket{ij}` denotes nearest neighbors on a square lattice in the xy-plane, :math:`t` are the hopping amplitudes, :math:`\theta_{ij}` are the Peierls phases, and :math:`c^{(\dagger)}` are the (creation)annihilation operators for spinless fermions :cite:`Andrews20`.

Working in Landau gauge in the x-direction, we can see that the Peierls phases are given as

.. math::
   \theta_{ij} = 2\pi n_\phi (Y_j - Y_i) \left( X_i + \frac{X_j - X_i}{2} \right),

where the derivation is given in the :doc:`Peierls phases <peierls_phases>` documentation. Hence, the complex phase factor acquired by each hopping describes a rectangular magnetic unit cell extended in the x-direction.

For example, for a flux density of :math:`n_\phi\equiv p/q=1/4`, the magnetic unit cell will be described by a rectangle of dimension 4x1, with constant x-hopping amplitudes of t, but varying y-hopping amplitudes of :math:`t`, :math:`t e^{\pi / 2}`,  :math:`t e^{\pi}`, :math:`t e^{2\pi 3/4}`, respectively.

Applying the Schrödinger equation to a wavefunction at a site indexed by (m,n) yields

.. math::
   E\Psi_{m,n} = t(\Psi_{m+1,n} + e^{i 2\pi n_\phi m}\Psi_{m, n+1}) + \mathrm{H.c.}.

Next, invoking a plane-wave ansatz for the wavefunction :math:`\Psi_{m,n}=e^{i k_x m} e^{i k_y n} \psi_m`, where :math:`0 \leq k_x < 2\pi/q` and :math:`0 \leq k_y < 2\pi`, yields

.. math::
   E\psi_{m} = t(e^{i k_x} \psi_{m+1} + e^{i (2\pi n_\phi m + k_y)}\psi_{m}) + \mathrm{H.c.}.

Hence, the equation may be written simply as

.. math::
   E\psi_{m} = B^*\psi_{m-1} + A_m \psi_{m} + B\psi_{m+1}

where

.. math::
   \begin{align}
       A_m &= 2\cos(2\pi n_\phi m + k_y), \\
       B &= t e^{i k_x}.
   \end{align}

Our wavefunction ansatz satisfies the generalized Bloch condition in the x-direction, given as :math:`\psi_{m+q}=\psi_{m}`.

For the case of :math:`n_\phi=1/4`, the :math:`q\times q` Hamiltonian matrix may be written as

.. math::
   \mathbf{H} =
   \begin{pmatrix}
   A_0 & B & 0 & B^* \\
   B^* & A_1 & B & 0 \\
   0 & B^* & A_2 & B \\
   B & 0 & B^* & A_3
   \end{pmatrix},

and hence there are 4 energy bands in the spectrum.

With the appropriate Peierls phases, this procedure may be readily generalized to other lattices, e.g. triangular, honeycomb, or kagome, and general flux densities :math:`n_\phi\equiv p/q`. We note that for the case of multiple sublattices, care needs to be taken for the generalized Bloch condition, e.g. for the honeycomb lattice with A and B sublattices, the generalized Bloch condition is :math:`\psi_{m+M}=e^{i k_x M b}\psi_{m}`, such that M=q(2q) for even(odd) p, and :math:`m\in[0,M]`.
