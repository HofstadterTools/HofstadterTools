Quantum Geometry
================

In this section, we will demonstrate how to numerically compute the quantum geometric tensor. :cite:`Parameswaran13, Claassen15`

The quantum geometric tensor is defined as

.. math::
       \mathcal{R}_{\alpha \beta}(\mathbf{k}) = \braket{\partial_{k_\alpha} u_\mathbf{k}|\partial_{k_\beta} u_\mathbf{k}} - \braket{\partial_{k_\alpha} u_\mathbf{k}| u_\mathbf{k}} \braket{u_\mathbf{k}| \partial_{k_\beta} u_\mathbf{k}},

where :math:`\ket{u_\mathbf{k}}` is the eigenvector at momentum :math:`\mathbf{k}`. The real part of the quantum geometric tensor is given by the Fubini-Study metric :math:`g_{ij}(\mathbf{k})=\Re[\mathcal{R}_{ij}]`, which corresponds to the distance between eigenstates on the Bloch sphere, whereas the imaginary part of the quantum geometric tensor is given by the Berry curvature :math:`\mathcal{B}(\mathbf{k})=-2 \Im [\mathcal{R}_{\alpha\beta}(\mathbf{k})]` with :math:`\alpha\neq\beta`.

As for the Chern numbers, there are several methods to compute the quantum geometric tensor. In analogy to the  :doc:`previous section <chern_numbers>`, we may start by naively writing down the numerical derivatives:

.. math::
   \bra{\partial_\alpha u_\mathbf{k}} = \lim_{\Delta_{\alpha}\to 0} \frac{\bra{ u_{\mathbf{k}+\Delta_\alpha}} - \bra{u_\mathbf{k}}}{\Delta_\alpha},

.. math::
   \ket{\partial_\beta u_\mathbf{k}} = \lim_{\Delta_{\beta}\to 0} \frac{\ket{ u_{\mathbf{k}+\Delta_\beta}} - \ket{u_\mathbf{k}}}{\Delta_\beta}.

Inserting these into the expression for the quantum geometric tensor and simplifying, yields

.. math::
   \mathcal{R}_{\alpha \beta}(\mathbf{k}) = \lim_{\Delta_\alpha\to 0} \lim_{\Delta_\beta \to 0} \frac{1}{\Delta_\alpha \Delta_\beta}\left[ \braket{u_{\mathbf{k}+\Delta_\alpha}| u_{\mathbf{k}+\Delta_\beta}} - \braket{u_{\mathbf{k}+\Delta_\alpha}| u_\mathbf{k}} \braket{ u_\mathbf{k}| u_{\mathbf{k}+\Delta_\beta}}\right].

This is a valid expression for the quantum geometric tensor, however it is gauge dependent and therefore care needs to be taken to ensure that the numerical derivatives are accurate. In practise, this often means choosing an interval for these derivatives that is many orders of magnitude smaller than the k mesh.

An alternative method for computing the quantum geometric tensor is to use projectors, such that

.. math::
   \mathcal{R}_{\alpha, \beta}(\mathbf{k}) = \partial_{k_\alpha}Q_\mathbf{k} \partial_{k_\beta} \mathcal{P}_\mathbf{k},

where :math:`\mathcal{P}_\mathbf{k} = \ket{u_\mathbf{k}} \bra{u_\mathbf{k}}` is the band projector and :math:`\mathcal{Q}_\mathbf{k}=1 - \mathcal{P}_\mathbf{k}` is the orthogonal band projector. This expression is now gauge invariant and therefore less susceptible to errors. The projector formalism can also be naturally extended to the case of band touching, as discussed in several works, e.g. Appendix C.4 of :cite:`Hirschmann23`.

Crucially, since band geometry and topology are components of the same tensor, we can derive relations between them, namely

.. math::
   \begin{align}
   \mathcal{D}(\mathbf{k})&=\text{det}\,g(\mathbf{k}) - \frac{1}{4}|\mathcal{B}(\mathbf{k})|^2 \geq 0, \\
   \mathcal{T}(\mathbf{k})&=\text{tr}\,g(\mathbf{k}) - |\mathcal{B}(\mathbf{k})| \geq 0,
   \end{align}

where we define :math:`\mathcal{D}` as the determinant inequality saturation measure (DISM) and :math:`\mathcal{T}` as the trace inequality saturation measure (TISM). It has been shown analytically that when the trace(determinant) inequality is saturated for a Chern band, the algebra of projected density operators is identical(isomorphic) to that in Landau levels :cite:`Roy14`. For the conventional Hofstadter model, this means that as we approach the continuum limit :math:`n_\phi\to 0`, the TISM and DISM will approach zero.
