Troubleshooting
===============

**Q.** I get an argparse error whenever one of my command line arguments is a negative number. What is correct syntax so that python does not interpret my negative numbers as flags?

**A.** If you would like to parse negative number arguments to command line flags, it's best to use quotation marks with a preceding space, e.g. ``-t 1 0 " -0.25"``.

**Q.** I have plotted a Hofstadter butterfly for some custom model but there are spurious straggeling bands and aperiodicity in the spectrum. What can I do to fix this?

**A.** By default, the flux density in HofstadterTools is defined with respect to the lattice unit cell area. However, in some models, the minimal plaquette around which a particle can hop encloses an area that is smaller than the unit cell area. In these cases, in order to both restore periodicity and view the complete butterfly spectrum, you may need to define the flux density with respect to the area of a minimal plaquette. In general, compute the ratio ``n`` of the effective unit cell area (spanned by the hopping terms) and the area of a minimal hopping plaquette, and then append the flag ``--periodicty n``.