Band Structures
===============

In this section, we will demonstrate how the band structures are computed for the conventional Hofstadter Model

.. math::
   H = -t \sum_{\braket{ij}} e^{i \theta_{ij}} c_i^\dagger c_j + \mathrm{H.c.},

where :math:`t` is the hopping amplitude, :math:`\braket{ij}` denotes nearest neighbors on a square lattice in the xy-plane with lattice constant :math:`a=1`, :math:`\theta_{ij}` are the Peierls phases acquired by hopping from site :math:`i=(X_i, Y_i)` to :math:`j=(X_j, Y_j)` due to a magnetic field in the z-direction, and :math:`c^{(\dagger)}` are the particle (creation)annihilation operators.

Working in Landau gauge in the y-direction, we have shown in the :doc:`previous section <peierls_phases>` that the Peierls phases may be written as

.. math::
   \theta_{ij} = -2\pi n_\phi \Delta X \left( Y_i + \frac{\Delta Y}{2} \right),

where :math:`\Delta X = X_j - X_i` and :math:`\Delta Y = Y_j - Y_i`. Hence, the complex phase factor acquired by each hopping describes a rectangular magnetic unit cell extended in the y-direction.

For example, for a rational flux density of :math:`n_\phi\equiv p/q=1/4`, the magnetic unit cell will be described by a rectangle of dimension :math:`1\times4`, with constant y-hopping amplitudes of :math:`-t`, but varying x-hopping amplitudes of :math:`-t`, :math:`-t e^{-\mathrm{i}\pi / 2}`,  :math:`-t e^{-\mathrm{i}\pi}`, :math:`-t e^{-\mathrm{i}3\pi /2}`, respectively.

Applying the Schrödinger equation to a wavefunction at a site indexed by :math:`(m,n)` yields

.. math::
   E\Psi_{m,n} = -t(e^{-\mathrm{i} 2\pi n_\phi n}\Psi_{m+1,n} + \Psi_{m, n+1}) + \mathrm{H.c.}.

Next, invoking a plane-wave ansatz for the wavefunction in the x-direction :math:`\Psi_{m,n}=e^{\mathrm{i}\mathbf{k}\cdot\mathbf{r}} \psi_n`, yields

.. math::
   E\psi_{n} = -t(e^{-\mathrm{i}2\pi n_\phi n + \mathrm{i}k_x}\psi_{n} + e^{\mathrm{i}k_y}\psi_{n+1}) + \mathrm{H.c.}.

Hence, the equation may be written simply as

.. math::
   E\psi_{n} = B^*\psi_{n-1} + A_n \psi_{n} + B\psi_{n+1}

where

.. math::
   \begin{align}
       A_m &= -2t\cos(2\pi n_\phi m - k_x), \\
       B &= -t e^{\mathrm{i} k_y}.
   \end{align}

This difference equation is known as the Harper equation and allows us to write the Schrödinger equation in matrix form, such that  :math:`H \boldsymbol{\psi} = E \boldsymbol{\psi}`, where :math:`H` is a :math:`q\times q` Hamiltonian matrix and :math:`\boldsymbol{\psi}=(\psi_0, \psi_1, \dots, \psi_{q-1})^\intercal` is a vector of length :math:`q`. :cite:`Harper55`

For the case of :math:`n_\phi=1/4`, the :math:`q\times q` Hamiltonian matrix may be written explicitly as

.. math::
   \mathbf{H} =
   \begin{pmatrix}
   A_0 & B & 0 & B^* \\
   B^* & A_1 & B & 0 \\
   0 & B^* & A_2 & B \\
   B & 0 & B^* & A_3
   \end{pmatrix},

and hence there are 4 energy bands in the spectrum.

With the appropriate Peierls phases, this procedure may be readily generalized to other lattices, e.g. triangular, honeycomb, or kagome, and general flux densities :math:`n_\phi\equiv p/q`.
