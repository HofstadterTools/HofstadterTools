Troubleshooting
===============

.. admonition:: Question

	I get an argparse error whenever one of my command line arguments is a negative number. What is correct syntax so that python does not interpret my negative numbers as flags?

If you would like to parse negative number arguments to command line flags, it's best to use quotation marks with a preceding space, e.g. ``-t 1 0 " -0.25"``.

.. admonition:: Question

	I have plotted a Hofstadter butterfly for some custom model but there are spurious straggeling bands and aperiodicity in the spectrum. What can I do to fix this?

By default, the flux density in HofstadterTools is defined with respect to the lattice unit cell area. However, in some models, the minimal plaquette around which a particle can hop encloses an area that is smaller than the unit cell area. In these cases, in order to both restore periodicity and view the complete butterfly spectrum, you may need to define the flux density with respect to the area of a minimal plaquette. In general, compute the ratio ``n`` of the effective unit cell area (spanned by the hopping terms) and the area of a minimal hopping plaquette, and then append the flag ``--periodicity n``.

.. admonition:: Question

	When computing for the ``square`` or ``triangular`` lattices, the flags ``-alpha`` and ``-theta`` do not seem to do anything. Why not?

The flags ``-alpha`` and ``-theta`` specify the anisotropy and obliqueness of the underlying Bravais lattice. In the case, of the ``square`` and ``triangular`` lattices, both of these values are fixed by definition and so cannot be changed. If you would like to change the anisotropy and obliqueness of a lattice with a single-site basis, you can use the ``bravais`` lattice parameter instead. For lattices with a multi-site basis, e.g. ``honeycomb`` or ``kagome``, you can use the ``-alpha`` and ``-theta`` flags as normal, since the basis sites are defined in fractional coordinates with respect to the Bravais vectors.

.. admonition:: Question

	I have run the code with the ``-save`` flag but I cannot find the output. Where is the output saved?

If you executed the code in the ``code`` directory (recommended), then the output will be saved by default in the ``data``, ``figs``, and ``logs`` directories. If you executed the code in any other location, then all of the output will be saved in the working directory.
