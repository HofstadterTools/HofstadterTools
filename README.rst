HofstadterTools
===============

|docs| |pytests| |license|

.. |docs| image:: https://github.com/HofstadterTools/HofstadterTools/actions/workflows/docs.yml/badge.svg
   :target: https://github.com/HofstadterTools/HofstadterTools/actions/workflows/docs.yml

.. |pytests| image:: https://github.com/HofstadterTools/HofstadterTools/actions/workflows/pytests.yml/badge.svg
   :target: https://github.com/HofstadterTools/HofstadterTools/actions/workflows/pytests.yml

.. |license| image:: https://badgen.net/badge/license/GPLv3/blue
   :target: https://www.gnu.org/licenses/gpl-3.0

* Paper: https://arxiv.org/abs/2311.18726
* Documentation: https://hofstadter.tools or `hof.tools <https://hofstadter.tools>`__
* GitHub Repository: https://github.com/HofstadterTools/HofstadterTools

**H**\ ofstadter\ **T**\ ools (\ **HT**) is a set of Python programs and classes for analyzing the Hofstadter model, which describes the behavior of non-interacting quantum particles hopping on a lattice coupled to a gauge field. This package can be used to compute the band structure of a *generalized* Hofstadter model on *any* regular Euclidean lattice, as well as its key properties, such as quantum geometry and topology.

.. admonition:: About the name

    Philip Harper first derived the difference equation for the model [`Harper55 <https://dx.doi.org/10.1088/0370-1298/68/10/304>`__], which was later analyzed in detail by Mark Azbel [`Azbel64 <http://jetp.ras.ru/cgi-bin/e/index/e/19/3/p634?a=list>`__], and finally plotted by Douglas Hofstadter [`Hofstadter76 <https://link.aps.org/doi/10.1103/PhysRevB.14.2239>`__]. Consequently, the formal name for the model is the *Harper-Azbel-Hofstadter model* to credit its three main contributors.

Quick Start
-----------

Using HofstadterTools is quick and easy!

.. code:: console

    pip install HofstadterTools

You can now access the ``band_structure`` and ``butterfly`` programs from any location.

.. code:: console

    band_structure -lat square -nphi 1 4 --display both --wilson

|image1| |image2| |image3|

.. |image1| image:: https://github.com/HofstadterTools/HofstadterTools/blob/2abdf3cb7c6ebfcce165e52a8020a329e5895313/docs/source/images/overview/band_structure_3D_both_square_nphi_1_4_t_1.png?raw=true
    :width: 32 %
    :alt: 3D Band Structure
.. |image2| image:: https://github.com/HofstadterTools/HofstadterTools/blob/2abdf3cb7c6ebfcce165e52a8020a329e5895313/docs/source/images/overview/wilson_both_square_nphi_1_4_t_1.png?raw=true
    :width: 32 %
    :alt: Wilson Loops
.. |image3| image:: https://github.com/HofstadterTools/HofstadterTools/blob/2abdf3cb7c6ebfcce165e52a8020a329e5895313/docs/source/images/overview/band_structure_2D_both_square_nphi_1_4_t_1.png?raw=true
    :width: 32 %
    :alt: 2D Band Structure

.. code:: console

    butterfly -lat square -q 97 --color point --wannier --plot_lattice

|image4| |image5| |image6|

.. |image4| image:: https://github.com/HofstadterTools/HofstadterTools/blob/2abdf3cb7c6ebfcce165e52a8020a329e5895313/docs/source/images/overview/butterfly_square_q_97_t_1_col_point_avron.png?raw=true
    :width: 32 %
    :alt: Butterfly Spectrum
.. |image5| image:: https://github.com/HofstadterTools/HofstadterTools/blob/2abdf3cb7c6ebfcce165e52a8020a329e5895313/docs/source/images/overview/wannier_square_q_97_t_1_col_point_avron.png?raw=true
    :width: 32 %
    :alt: Wannier Diagram
.. |image6| image:: https://github.com/HofstadterTools/HofstadterTools/blob/2abdf3cb7c6ebfcce165e52a8020a329e5895313/docs/source/images/overview/lattice.png?raw=true
    :width: 32 %
    :alt: Lattice

Voilà! You have just plotted the Hofstadter band structure for nearest-neighbor hopping on the square lattice at flux density :math:`n_\phi=1/4`, together with the corresponding butterfly spectrum at :math:`q=97`. You can append ``--help`` to either of these programs to view the list of options. Alternatively, you can explore the `gallery <https://hofstadter.tools/gallery.html>`__ and `code reference <https://hofstadter.tools/_autosummary/functions.html>`__ to see what HofstadterTools has to offer.

Installation
------------

This package was developed using Ubuntu 20.04.6 (x86_64) with Python=3.10.13, however it is designed to be platform-independent and can work with any Python>=3.9.

Basic install
~~~~~~~~~~~~~

For basic usage of HofstadterTools, i.e. in cases where you *do not* plan on editing the source code, you can install the package quickly and easily from PyPI.

.. code:: console

    pip install HofstadterTools

To avoid dependency clashes, we recommend the use of a python virtual environment, such as ``venv``.

.. warning::

    If you install HofstadterTools into a conda environment, you may see a ``libGL error`` when you run the programs. This is a known problem with the ``libstdc++.so`` file in Conda and should not affect the functionality of HofstadterTools.

Advanced install
~~~~~~~~~~~~~~~~

For advanced usage of HofstadterTools, i.e. in cases where you *do* plan on editing the source code, you can install the package from source.

.. code:: console

    user@domain:any/path$ git clone git@github.com:HofstadterTools/HofstadterTools.git
    user@domain:any/path$ cd HofstadterTools
    user@domain:path/to/HofstadterTools$ pip install -e .

The optional ``-e`` flag indicates an editable install. Alternatively, if you plan on building the documentation locally, the optional ``docs`` dependencies need to be installed.

.. code:: console

    user@domain:path/to/HofstadterTools$ pip install -e ".[docs]"

.. note::

    Implementing custom lattices with more than one site per unit cell requires an advanced install.

Testing
~~~~~~~

You can confirm that HofstadterTools is correctly installed by running the pytests.

.. code:: console

    user@domain:any/path$ pytest --pyargs HT

Once the *project* ``HofstadterTools`` is installed, the *package* ``HT`` will be available in your python environment. In addition, you can access the programs ``band_structure``, ``butterfly``, ``plot_band_structure``, and ``plot_butterfly``, from any location.

.. code:: console

    user@domain:any/path$ band_structure --help
    user@domain:any/path$ butterfly --help
    user@domain:any/path$ plot_band_structure --help
    user@domain:any/path$ plot_butterfly --help

The ``plot_*`` programs are used to replot band_structures / butterflies that have been saved to file.

Directory Structure
-------------------

* **src** -- sources root with the ``HT`` package, along with its configuration settings, subpackages, and programs. A detailed description of the available `programs <https://hofstadter.tools/tutorials.html>`__ and `namespace packages <https://hofstadter.tools/_autosummary/functions.html>`__ is in the documentation.

  * **HT** -- ``HT`` package

    * **configuration** -- user-defined configuration files for the programs.
    * **functions** -- helper functions for the programs.
    * **models** -- model classes for the programs.
    * **plot** -- location of the plot scripts.
    * **tests** -- unit tests for the programs.

* **data** -- output destination for raw data files (if programs are run explicitly from their file location, otherwise the output destination is the current working directory).

  * **band_structure** -- data generated by the band_structure program.
  * **butterfly** -- data generated by the butterfly program.

* **docs** -- location of the sphinx documentation. To view the documentation locally, compile by running ``make html`` or ``make clean html`` and then open ``build/html/index.html`` in a web browser. This assumes that the optional ``docs`` dependencies are installed.

  * **build** -- compiled documentation (once built).
  * **source** -- documentation source.

* **figs** -- output destination for the figures (if programs are run explicitly from their file location, otherwise the output destination is the current working directory).

  * **band_structure** -- figures generated by the band_structure program.
  * **butterfly** -- figures generated by the butterfly program.

* **logs** -- output destination for the log files (if programs are run explicitly from their file location, otherwise the output destination is the current working directory).

  * **band_structure** -- logs generated by the band_structure program.
  * **butterfly** -- logs generated by the butterfly program.

* **paper** -- summary paper introducing HofstadterTools. The formatted pdf can be downloaded as an ``artifact`` of the ``production-pdf`` workflow under the GitHub actions tab.

How to Cite
-----------

If you have found HofstadterTools useful, it would be greatly appreciated if you could cite us in your work. Please find the bibtex reference below.

.. code-block:: bibtex

  @misc{HofstadterTools,
  title={HofstadterTools: A Python package for analyzing the Hofstadter model},
  author={Bartholomew Andrews},
  year={2023},
  eprint={2311.18726},
  archivePrefix={arXiv},
  primaryClass={cond-mat.mes-hall}
  }

Acknowledgments
---------------

We thank Gunnar Möller, Titus Neupert, Rahul Roy, Alexey Soluyanov, Michael Zaletel, Daniel Parker, Stefan Divic, Johannes Mitscherling, and Mathi Raja, for useful discussions. This project was funded by the Swiss National Science Foundation under Grant No. `P500PT_203168 <https://data.snf.ch/grants/grant/203168>`__, and supported by the U.S. Department of Energy, Office of Science, Basic Energy Sciences, under Early Career Award No. DE-SC0022716.

Contributing
------------

The Hofstadter model is an active field of research and therefore HofstadterTools will never be complete. Here is a list of some features that we have on the pipeline to be implemented (in no particular order):

* support for hyperbolic lattices [`Stegmaier22 <https://link.aps.org/doi/10.1103/PhysRevLett.128.166402>`__]
* support for fractal lattices [`Chen20 <https://doi.org/10.1007/s00220-020-03850-w>`__]
* support for higher-dimensional lattices [`DiColandrea22 <https://dx.doi.org/10.1088/1367-2630/ac4126>`__]
* support for quasicrystals [`Ghadimi22 <https://link.aps.org/doi/10.1103/PhysRevB.106.L201113>`__]
* support for open boundary conditions [`Pena23 <https://doi.org/10.1016/j.rinp.2023.106257>`__]
* interface to quantum chemistry codes [`Bodesheim23 <https://doi.org/10.1038/s41699-023-00378-0>`__]
* capability to compute the non-Abelian `Hofstadter moth` [`Osterloh05 <https://link.aps.org/doi/10.1103/PhysRevLett.95.010403>`__], [`Yang20 <https://doi.org/10.1038/s41377-020-00384-7>`__]
* capability to compute Chern numbers using bulk-edge correspondence [`Agazzi14 <https://doi.org/10.1007/s10955-014-0992-0>`__]
* capability to generate the potential function corresponding to hopping amplitudes [`Yilmaz17 <https://link.aps.org/doi/10.1103/PhysRevA.95.063628>`__]
* implementation of other topological flat-band models for benchmarking (e.g. chiral pi-flux model) [`Neupert11 <https://link.aps.org/doi/10.1103/PhysRevLett.106.236804>`__]

Contributions are always welcome! The easiest way to contribute is to submit a pull request on `GitHub <https://github.com/HofstadterTools/HofstadterTools>`__ or contact `Bart Andrews <https://bartandrews.me>`__ if you have any feedback.
