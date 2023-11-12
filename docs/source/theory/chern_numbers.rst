Chern Numbers
=============

In this section, we will demonstrate how to numerically compute the first Chern number. :cite:`AidelsburgerPhD, SoluyanovPhD, Goudarzi22`

The Chern number is defined as

.. math::
   C = \frac{1}{2\pi} \int_\mathrm{BZ} \mathcal{B}(\mathbf{k}) \mathrm{d}^2 k,

where :math:`\mathcal{B}(\mathbf{k})=\epsilon^{\alpha \beta}\partial_{k_\alpha} \mathbf{A}_{\beta}(\mathbf{k})` with :math:`\alpha\neq\beta` is the Berry curvature, :math:`\mathbf{A}_\beta(\mathbf{k})=-\mathrm{i}\braket{u_\mathbf{k}|\partial_{k_\beta}u_\mathbf{k}}` is the Berry connection, and the integral is taken over the Brillouin zone.

Writing the Berry curvature out in full yields

.. math::
   \mathcal{B}(\mathbf{k}) = -\mathrm{i} (\partial_{k_\alpha} \braket{u_\mathbf{k}|\partial_{k_\beta}u_\mathbf{k}} - \partial_{k_\beta} \braket{u_\mathbf{k}|\partial_{k_\alpha}u_\mathbf{k}}).

Using the fact that the numerical derivative may be written as

.. math::
   \partial_{k_\alpha} u_\mathbf{k} = \lim_{\Delta_\alpha\to 0} \frac{u_{\mathbf{k}+\Delta_\alpha}-u_\mathbf{k}}{\Delta_\alpha},

where :math:`\Delta_\alpha=\Delta k_\alpha`, we may express the double derivative as

.. math::
   \partial_{k_\alpha} \braket{u_\mathbf{k}|\partial_{k_\beta} u_\mathbf{k}} = \lim_{\Delta_\alpha\to 0}\lim_{\Delta_\beta\to 0} \frac{\braket{u_{\mathbf{k}+\Delta_\alpha}| u_{\mathbf{k}+\Delta_\beta+\Delta_\alpha}}-\braket{u_\mathbf{k}|u_{\mathbf{k}+\Delta_\beta}}}{\Delta_\alpha \Delta_\beta}.

Inserting this into the expression for the Berry curvature yields

.. math::
   \mathcal{B}(\mathbf{k}) = -\mathrm{i} \lim_{\Delta_\alpha\to 0}\lim_{\Delta_\beta\to 0} \frac{1}{\Delta_\alpha\Delta_\beta} \left( \braket{u_{\mathbf{k}+\Delta_\alpha}| u_{\mathbf{k}+\Delta_\beta+\Delta_\alpha}}-\braket{u_\mathbf{k}|u_{\mathbf{k}+\Delta_\beta}} - \braket{u_{\mathbf{k}+\Delta_\beta}| u_{\mathbf{k}+\Delta_\alpha+\Delta_\beta}}+\braket{u_\mathbf{k}|u_{\mathbf{k}+\Delta_\alpha}} \right).


However, it is difficult to compute the Berry curvature numerically using this formula because it is not gauge invariant. Instead, exploiting the property that the Berry curvature is real, we may write it as a product of Berry phases around a plaquette

.. math::
   \mathcal{B}(\mathbf{k}) = - \lim_{\Delta_\alpha\to 0}\lim_{\Delta_\beta\to 0} \frac{1}{\Delta_\alpha \Delta_\beta} \text{Im}\;\log (\braket{u_{\mathbf{k}}| u_{\mathbf{k}+\Delta_\alpha}}\braket{u_{\mathbf{k}+\Delta_\alpha}| u_{\mathbf{k}+\Delta_\alpha+\Delta_\beta}}\braket{u_{\mathbf{k}+\Delta_\alpha+\Delta_\beta}| u_{\mathbf{k}+\Delta_\beta}} \braket{u_{\mathbf{k}+\Delta_\beta}| u_{\mathbf{k}}}).

This expression is gauge invariant and the Chern number is then computed by numerically integrating the Berry curvature at each plaquette and dividing by :math:`2\pi`. In the limit of small plaquettes, we recover an integer. This equation is also known as the Fukui formula. :cite:`Fukui05`

In the case of band touching, each phase factor is replaced by the determinant of the phase factor matrix. For example, in the case of two bands touching, the phase factors

.. math::
   \braket{u_{\mathbf{k}}| u_{\mathbf{k}+\Delta_\alpha}} \to
   \begin{Vmatrix}
   \braket{u_{0, \mathbf{k}}| u_{0, \mathbf{k}+\Delta_\alpha}} & \braket{u_{0, \mathbf{k}}| u_{1, \mathbf{k}+\Delta_\alpha}} \\
   \braket{u_{1, \mathbf{k}}| u_{0, \mathbf{k}+\Delta_\alpha}} & \braket{u_{1, \mathbf{k}}| u_{1, \mathbf{k}+\Delta_\alpha}}
   \end{Vmatrix},

etc., may be used to compute the Chern number of the combined band pair. The Chern numbers in a Hofstadter spectrum always sum to zero.
