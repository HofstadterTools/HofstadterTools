Troubleshooting
===============

General Questions
-----------------

.. admonition:: Question

	I get an argparse error whenever one of my command line arguments is a negative number. What is correct syntax so that python does not interpret my negative numbers as flags?

If you would like to parse negative number arguments to command line flags, it's best to use quotation marks with a preceding space, e.g. ``-t 1 0 " -0.25"``.

.. admonition:: Question

	When computing for the ``square`` or ``triangular`` lattices, the flags ``-alpha`` and ``-theta`` do not seem to do anything. Why not?

The flags ``-alpha`` and ``-theta`` specify the anisotropy and obliqueness of the underlying Bravais lattice. In the case of the ``square`` and ``triangular`` lattices, both of these values are fixed by definition and so cannot be changed. If you would like to change the anisotropy and obliqueness of a lattice with a single-site basis, you can use the ``bravais`` lattice parameter instead. For lattices with a multi-site basis, e.g. ``honeycomb`` or ``kagome``, you can use the ``-alpha`` and ``-theta`` flags as normal, since the basis sites are defined in fractional coordinates with respect to the Bravais vectors.

.. admonition:: Question

	I have run the code with the ``-save``/``-log`` flag but I cannot find the output. Where is the output saved?

If you executed the code in the ``code`` directory (recommended), then the output will be saved by default in the ``data``, ``figs``, and ``logs`` directories. If you executed the code in any other location, then all of the output will be saved in the current working directory.

.. admonition:: Question

	I would like to implement a custom lattice, which has a multi-site basis that is different than ``honeycomb`` or ``kagome``. How can I do this?

In HofstadterTools, the ``custom`` lattice parameter is used for this purpose. Navigate to the ``unit_cell`` method of the ``Hofstadter`` class in ``models/hofstadter.py`` and enter your desired set of basis vectors and high-symmetry points under the ``elif self.lat == "custom":`` clause. By default, the ``custom`` lattice is identical to ``kagome``.

.. admonition:: Question

	I would like to input a list of hopping parameters from a file. How can I do this?

If you are executing your code from within the ``code`` directory (recommended), then you can edit the configuration file ``hopping_input.txt`` located in ``code/configuration``, and then append the flag ``-input`` to your program. The file consists of rows with two columns, separated by a tab. The left column is the NN group and the right column is the hopping amplitude. For example, in order to input the flags ``-t 1 0 " -0.25"``, we can create the following file:

.. code:: console

		1 1
		3 -0.25

If you are executing your program from a different location, then the ``-input`` flag will search for a file called ``hopping_input.txt`` in the current working directory.

.. admonition:: Question

	The progress bar becomes much shorter when I run using the ``-log`` flag. Why is this?

By default, the tqdm progress bar is streamed to stderr. Hence, the progress bar is shortened to 10 characters when logging to make the log files more readable.

.. admonition:: Question

	I have adjusted the ``-alpha`` lattice anisotropy flag and now the band structures and butterfly spectra look odd. Is this expected behavior?

In HofstadterTools, nth-nearest neighbors are defined as a group of sites that are the same radius from a reference site. We stress that these points lie on a *circle* and not, for example, on an ellipse that is scaled by the Bravais lattice vector lengths. Hence, by adjusting the lattice anisotropy, you may inadvertently change the number of nearest neighbors.

band_structure Questions
------------------------

.. admonition:: Question

	I would like to analyze the quantum geometry of the bands, but I do not see any quantum geometry data in the output table of the ``band_structure`` program. How can I do this?

By default, the quantum geometry computations are turned off in the ``band_structure`` program (for speed reasons). If you would like to turn these on, or configure the output table in any other way, you can edit the table column selector at ``code/configuration/band_structure.py``. Note that the table columns are grouped by computational expense. For example, if you would like to output the normalized Fubini-Study metric fluctuations ``std_g_norm``, then you could also output the Brillouin-zone-averaged trace inequality saturation measure ``T``  at negligible additional cost.

.. admonition:: Question

	I would like to compute the band structure of my custom Hamiltonian but it is running very slowly compared to standard examples. Why is this?

The code for constructing a Hamiltonian matrix for a generalized Hofstadter model on any regular Euclidean lattice is expensive. In light of this, we have hardcoded the most common Hofstadter Hamiltonians, i.e. the Hofstadter Hamiltonians on the square/triangular/honeycomb/kagome lattices with nearest-neighbor interactions. In all other cases, the generic Hamiltonian constructor will be called. If you are interested in one custom Hamiltonian in particular, and really need to compute its complete band structure more quickly, then consider adding it to ``functions/models.py`` in a similar format to the other hardcoded Hamiltonians, e.g. ``BasicKagomeHamiltonian``. We note that the ``butterfly`` program does not suffer from this issue, since the diagonalization is performed only at a single :math:`k` point.

.. admonition:: Question

	I have computed the band structure for a particular model and I have noticed that certain bands are not "touching" when they should be, or visa versa. How can I fix this?

Due to the discrete nature of the :math:`k` mesh, it is difficult to declare that certain bands are touching. For this purpose, HofstadterTools uses the band gap threshold flag ``-bgt``, which declares bands as touching when they are within this value of each other. If you notice that certain bands should/should not be touching, e.g. by noticing that the Chern numbers do not sum to zero, or you are simply suspicious of bands that are in close proximity, you can try decreasing the mesh size using the ``-samp`` flag and tweaking this ``-bgt`` value.

.. admonition:: Question

	I have computed the band structure for a kagome/custom lattice and it looks incorrect. Why is this?

When computing the complete band structure, it may be more difficult to spot when the ``--periodicity`` flag needs to be set. If in doubt, compute the corresponding butterfly spectrum and make sure that it has the correct periodicity.

butterfly Questions
-------------------

.. admonition:: Question

	I have plotted a Hofstadter butterfly for some custom model but there are spurious straggeling bands and aperiodicity in the spectrum. What can I do to fix this?

By default, the flux density in HofstadterTools is defined with respect to the lattice unit cell area. However, in some models, the minimal plaquette around which a particle can hop encloses an area that is smaller than the unit cell area. In these cases, in order to both restore periodicity and view the complete butterfly spectrum, you may need to define the flux density with respect to the area of a minimal plaquette. In general, compute the ratio ``n`` of the effective unit cell area (spanned by the hopping terms) and the area of a minimal hopping plaquette, and then append the flag ``--periodicity n``.

.. admonition:: Question

	I have plotted a Hofstadter butterfly using the ``--color`` flag and the code runs surprisingly quickly. How are the Chern numbers computed?

All of the Hofstadter butterflies are colored using the Streda-Widom Diophantine relation (see Appendix C of :cite:`DiColandrea22` for a derivation). We note that although the formula can unambiguously determine the Chern numbers for the case of the rectangular lattice, the natural window condition is not uniquely resolved in general. This *may* lead to minor imperfections in the coloring for other lattices, especially when plotting with an extremely high resolution, as scrutinized in Fig.4 of :cite:`Agazzi14` or :cite:`Avron14`. At the time of writing, there is no Diophantine equation that can uniquely determine the Chern numbers in the general case. For the ``--color`` flag, we make the choice of sacrificing precision for the sake of efficiency.

.. admonition:: Question

	I am trying to plot a plane-colored Hofstadter butterfly with high resolution but I find strange interpolated blobs in the fine structure of the spectrum. How can I fix this?

This is an indication that the dpi of the image is too low. Assuming that you have saved the output data for such a high-resolution spectrum (recommended), you can overwrite the ``args['dpi']`` parameter in the ``plot/butterfly.py`` script and try plotting it again. By default, the dpi is set to 300. This works reasonably well for :math:`M` values up to about 300, where :math:`M` is the number of bands in the spectrum. In general, we recommend setting a dpi value of greater than :math:`M` for best results.