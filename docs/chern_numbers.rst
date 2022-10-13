Chern numbers
=============

In this section, we will demonstrate how to numerically compute the Chern number.

The Chern number is defined as

.. math::

   C = \frac{1}{2\pi} \int_\mathrm{BZ} \mathcal{B}(\mathbf{k}) \mathrm{d}^2 k,

where :math:`\mathcal{B}(\mathbf{k})=\epsilon^{\alpha \beta}\partial_{k_\alpha} \mathbf{A}_{\beta}(\mathbf{k})` is the Berry curvature, and :math:`\mathbf{A}_\beta(\mathbf{k})=-i<u_\mathbf{k}|\partial_{k_\beta}u_\mathbf{k}>` is the Berry connection.

Writing the Berry connection out in full yields

.. math::

   \mathcal{B}(\mathbf{k}) = -i (\partial_{k_\alpha} <u_\mathbf{k}|\partial_{k_\beta}u_\mathbf{k}> - \partial_{k_\beta} <u_\mathbf{k}|\partial_{k_\alpha}u_\mathbf{k}>).

Using the fact that the numerical derivative may be written as

.. math::

   \partial_{k_\alpha} u_\mathbf{k} = \lim_{\Delta_\alpha\to 0} \frac{u_{\mathbf{k}+\Delta_\alpha}-u_\mathbf{k}}{\Delta_\alpha},

where :math:`\Delta_\alpha=\Delta k_\alpha`, we may express the double derivative as

.. math::

   \partial_{k_\alpha} <u_\mathbf{k}|\partial_{k_\beta} u_\mathbf{k}> = \lim_{\Delta_\alpha\to 0}\lim_{\Delta_\beta\to 0} \frac{<u_{\mathbf{k}+\Delta_\alpha}| u_{\mathbf{k}+\Delta_\beta+\Delta_\alpha}>-<u_\mathbf{k}|u_{\mathbf{k}+\Delta_\beta}>}{\Delta_\alpha \Delta_\beta}.

Inserting this into the expression for the Berry curvature yields

.. math::

   \mathcal{B}(\mathbf{k}) = -i \lim_{\Delta_\alpha\to 0}\lim_{\Delta_\beta\to 0} \frac{1}{\Delta_\alpha\Delta_\beta} \left( <u_{\mathbf{k}+\Delta_\alpha}| u_{\mathbf{k}+\Delta_\beta+\Delta_\alpha}>-<u_\mathbf{k}|u_{\mathbf{k}+\Delta_\beta}> - <u_{\mathbf{k}+\Delta_\beta}| u_{\mathbf{k}+\Delta_\alpha+\Delta_\beta}>+<u_\mathbf{k}|u_{\mathbf{k}+\Delta_\alpha}> \right).

Moreover, exploiting the property that the Berry curvature is real, we may write it as a sum of Berry phases around a plaquette

.. math::

   \mathcal{B}(\mathbf{k}) = - \lim_{\Delta_\alpha\to 0}\lim_{\Delta_\beta\to 0} \frac{1}{\Delta_\alpha \Delta_\beta} \text{Im}\;\log <u_{\mathbf{k}}| u_{\mathbf{k}+\Delta_\alpha}><u_{\mathbf{k}+\Delta_\alpha}| u_{\mathbf{k}+\Delta_\alpha+\Delta_\beta}><u_{\mathbf{k}+\Delta_\alpha+\Delta_\beta}| u_{\mathbf{k}+\Delta_\beta}> <u_{\mathbf{k}+\Delta_\beta}| u_{\mathbf{k}}>.

The Chern number is then computed by numerically integrating the Berry curvature at each plaquette and dividing by :math:`2\pi`. In the limit of small plaquettes, we recover an integer.

In the case of band touching, each phase factor is replaced by the determinant of the phase factor matrix. For example, in the case of two bands touching, the phase factors

.. math::

   <u_{\mathbf{k}}| u_{\mathbf{k}+\Delta_\alpha}> \to
   \begin{vmatrix}
   <u_{0, \mathbf{k}}| u_{0, \mathbf{k}+\Delta_\alpha}> & <u_{0, \mathbf{k}}| u_{1, \mathbf{k}+\Delta_\alpha}> \\
   <u_{1, \mathbf{k}}| u_{0, \mathbf{k}+\Delta_\alpha}> & <u_{1, \mathbf{k}}| u_{1, \mathbf{k}+\Delta_\alpha}>
   \end{vmatrix},

etc., may be used to compute the Chern number of the combined band pair.
