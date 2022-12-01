Fubini-Study Metric
===================

In this section, we will demonstrate how to numerically compute the Fubini-Study metric :cite:`Parameswaran13, Claassen15`.

The Fubini-Study metric is defined as

.. math::
   g_{\alpha \beta}(\mathbf{k}) = \frac{1}{2} \left[ \braket{\partial_{k_\alpha} u_\mathbf{k}|\partial_{k_\beta} u_\mathbf{k}} - \braket{\partial_{k_\alpha} u_\mathbf{k}| u_\mathbf{k}} \braket{u_\mathbf{k}| \partial_{k_\beta} u_\mathbf{k}} + (\alpha \leftrightarrow \beta)\right].

In analogy to the computation of the Chern number, outlined in the :doc:`Chern numbers <chern_numbers>` documentation, we may write down the numerical derivatives:

.. math::
   \bra{\partial_\alpha u_\mathbf{k}} = \lim_{\Delta_{\alpha}\to 0} \frac{\bra{ u_{\mathbf{k}+\Delta_\alpha}} - \bra{u_\mathbf{k}}}{\Delta_\alpha},

.. math::
   \ket{\partial_\beta u_\mathbf{k}} = \lim_{\Delta_{\beta}\to 0} \frac{\ket{ u_{\mathbf{k}+\Delta_\beta}} - \ket{u_\mathbf{k}}}{\Delta_\beta}.

Inserting these into the expression for the Fubini-Study metric and simplifying, yields

.. math::
   g_{\alpha \beta}(\mathbf{k}) = \frac{1}{2} \lim_{\Delta_\alpha\to 0} \lim_{\Delta_\beta \to 0} \frac{1}{\Delta_\alpha \Delta_\beta}\left[ \braket{u_{\mathbf{k}+\Delta_\alpha}| u_{\mathbf{k}+\Delta_\beta}} - \braket{u_{\mathbf{k}+\Delta_\alpha}| u_\mathbf{k}} \braket{ u_\mathbf{k}| u_{\mathbf{k}+\Delta_\beta}} + (\alpha \leftrightarrow \beta)\right].
