Fubini-Study Metric
===================

In this section, we will demonstrate how to numerically compute the Fubini-Study metric.

The Fubini-Study metric is defined as

.. math::

   g_{\alpha \beta}(\mathbf{k}) = \frac{1}{2} \left[ <\partial_{k_\alpha} u_\mathbf{k}|\partial_{k_\beta} u_\mathbf{k}> - <\partial_{k_\alpha} u_\mathbf{k}| u_\mathbf{k}>< u_\mathbf{k}| \partial_{k_\beta} u_\mathbf{k}> + (\alpha \leftrightarrow \beta)\right].

In analogy to the computation of the Chern number, outlined in the ``chern_numbers`` documentation, we may write down the numerical derivatives:

.. math::

   < \partial_\alpha u_\mathbf{k} | = \lim_{\Delta_{\alpha}\to 0} \frac{< u_{\mathbf{k}+\Delta_\alpha}| - <u_\mathbf{k}|}{\Delta_\alpha},

.. math::

   | \partial_\beta u_\mathbf{k} > = \lim_{\Delta_{\beta}\to 0} \frac{| u_{\mathbf{k}+\Delta_\beta}> - |u_\mathbf{k}>}{\Delta_\beta}.

Inserting these into the expression for the Fubini-Study metric and simplifying, yields

.. math::

   g_{\alpha \beta}(\mathbf{k}) = \frac{1}{2} \lim_{\Delta_\alpha\to 0} \lim_{\Delta_\beta \to 0} \frac{1}{\Delta_\alpha \Delta_\beta}\left[ <u_{\mathbf{k}+\Delta_\alpha}| u_{\mathbf{k}+\Delta_\beta}> - <u_{\mathbf{k}+\Delta_\alpha}| u_\mathbf{k}>< u_\mathbf{k}| u_{\mathbf{k}+\Delta_\beta}> + (\alpha \leftrightarrow \beta)\right].
