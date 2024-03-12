HofstadterTools
===============

|docs| |pytests| |pypi| |conda-forge| |license|

.. |docs| image:: https://github.com/HofstadterTools/HofstadterTools/actions/workflows/docs.yml/badge.svg
   :target: https://github.com/HofstadterTools/HofstadterTools/actions/workflows/docs.yml

.. |pytests| image:: https://github.com/HofstadterTools/HofstadterTools/actions/workflows/pytests.yml/badge.svg
   :target: https://github.com/HofstadterTools/HofstadterTools/actions/workflows/pytests.yml

.. |pypi| image:: https://img.shields.io/pypi/v/HofstadterTools
   :target: https://pypi.org/project/HofstadterTools/

.. |conda-forge| image:: https://img.shields.io/conda/v/conda-forge/hofstadtertools?label=conda-forge
   :target: https://anaconda.org/conda-forge/hofstadtertools

.. |license| image:: https://img.shields.io/github/license/HofstadterTools/HofstadterTools
   :target: https://www.gnu.org/licenses/gpl-3.0

* Paper: https://arxiv.org/abs/2311.18726
* Documentation: https://hofstadter.tools or `hof.tools <https://hofstadter.tools>`__
* GitHub Repository: https://github.com/HofstadterTools/HofstadterTools

**H**\ ofstadter\ **T**\ ools (\ **HT**) is a set of Python programs and classes for analyzing the Hofstadter model, which describes the behavior of non-interacting quantum particles hopping on a lattice coupled to a gauge field. This package can be used to compute the band structure of a *generalized* Hofstadter model on *any* regular Euclidean lattice, as well as its key properties, such as quantum geometry and topology.

.. admonition:: About the name

    Philip Harper first derived the difference equation for the model [`Harper55 <https://dx.doi.org/10.1088/0370-1298/68/10/304>`__], which was later analyzed in detail by Mark Azbel [`Azbel64 <http://jetp.ras.ru/cgi-bin/e/index/e/19/3/p634?a=list>`__], and finally plotted by Douglas Hofstadter [`Hofstadter76 <https://link.aps.org/doi/10.1103/PhysRevB.14.2239>`__]. Consequently, the formal name for the model is the *Harper-Azbel-Hofstadter model* to credit its three main contributors.

Quick Start
-----------

Installing HofstadterTools is quick and easy!

Using pip:

.. code:: console

    pip install HofstadterTools

Using conda:

.. code:: console

    conda install conda-forge::HofstadterTools

You can now access the ``band_structure`` and ``butterfly`` programs from any location. The ``band_structure`` program evaluates the Hofstadter band structure at fixed flux density :math:`n_\phi=p/q` for a range of momenta :math:`\mathbf{k}=(k_x,k_y)`, whereas the ``butterfly`` program evaluates the energy spectrum at fixed momentum :math:`\mathbf{k}=\mathbf{0}` for a range of :math:`n_\phi=p/q`, where :math:`p` and :math:`q` are coprime integers.

.. code:: console

    band_structure -lat square -nphi 1 4 --display both --wilson

|image1| |image2| |image3|

.. |image1| image:: https://github.com/HofstadterTools/HofstadterTools/blob/2abdf3cb7c6ebfcce165e52a8020a329e5895313/docs/source/images/overview/band_structure_3D_both_square_nphi_1_4_t_1.png?raw=true
    :width: 32 %
    :alt: 3D Band Structure
.. |image2| image:: https://github.com/HofstadterTools/HofstadterTools/blob/2abdf3cb7c6ebfcce165e52a8020a329e5895313/docs/source/images/overview/wilson_both_square_nphi_1_4_t_1.png?raw=true
    :width: 32 %
    :alt: Wilson Loops
.. |image3| image:: https://github.com/HofstadterTools/HofstadterTools/blob/fb6764269db9bfb84cf2f0fc7be0a729799db1bc/docs/source/images/overview/band_structure_2D_both_square_nphi_1_4_t_1.png?raw=true
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

For basic usage of HofstadterTools, i.e. in cases where you *do not* plan on editing the source code, you can install the package from a distribution.

Using pip
*********

1. [Optional] Create and activate a new venv environment. In the example below, ``my_env_name`` is the name of the venv and ``my_env_folder`` is the name of its folder. If needed, you can replace ``python3`` with ``python3.xx`` below to create a venv pinned to a particular Python version.

.. code:: console

    user@domain:any/path$ python3 -m my_env_name my_env_folder
    user@domain:any/path$ source path/to/my_env_folder/bin/activate

2. Install HofstadterTools from PyPI.

.. code:: console

    (my_env_name) user@domain:any/path$ pip install HofstadterTools

3. [Optional] Upgrade an existing HofstadterTools installation.

.. code:: console

    (my_env_name) user@domain:any/path$ pip install --upgrade HofstadterTools

You can verify the installation by typing ``pip list | grep HofstadterTools``, you can uninstall by typing ``pip uninstall HofstadterTools``, and you can deactivate the environment by typing ``deactivate``. The entire environment can be removed by deleting ``my_env_folder``.

Using conda
***********

1. [Optional] Create and activate a new conda environment. In the example below, ``my_env_name`` is the name of the conda environment. If needed, you can replace ``python=3`` with ``python=3.xx`` below to create a conda environment with a particular Python version pre-installed.

.. code:: console

    user@domain:any/path$ conda create -n my_env_name python=3
    user@domain:any/path$ conda activate my_env_name

2. Install HofstadterTools from conda-forge.

.. code:: console

    (my_env_name) user@domain:any/path$ conda install conda-forge::HofstadterTools

3. [Optional] Update an existing HofstadterTools installation.

.. code:: console

    (my_env_name) user@domain:any/path$ conda update HofstadterTools

You can verify the installation by typing ``conda list | grep hofstadtertools``, you can uninstall by typing ``conda remove HofstadterTools``, and you can deactivate the environment by typing ``conda deactivate``. The entire environment can be removed by typing ``conda remove -n my_env_name --all``.

.. warning::

    If you pip install HofstadterTools into a conda environment, you may see a ``libGL error`` when you run the programs. This is a known problem with the ``libstdc++.so`` file in Conda and should not affect the functionality of HofstadterTools.

Advanced install
~~~~~~~~~~~~~~~~

For advanced usage of HofstadterTools, i.e. in cases where you *do* plan on editing the source code, you can install the package from source.

1. Clone the HofstadterTools repository.

.. code:: console

    user@domain:any/path$ git clone git@github.com:HofstadterTools/HofstadterTools.git

2. Using pip, install the HofstadterTools package. This step can also be done in a virtual environment. The optional ``-e`` flag below indicates an editable install.

.. code:: console

    user@domain:path/to/HofstadterTools$ pip install -e .

Alternatively, if you plan on building the documentation locally, the optional ``docs`` dependencies need to be installed.

.. code:: console

    user@domain:path/to/HofstadterTools$ pip install -e ".[docs]"

3. [Optional] Build and view the documentation locally. The optional ``clean`` argument below removes files from the build directory, and ``firefox`` can be replaced with any web browser.

.. code:: console

    user@domain:path/to/HofstadterTools/docs$ make clean html
    user@domain:path/to/HofstadterTools/docs$ firefox build/html/index.html &

.. note::

    Building the documentation locally with the inheritance diagrams requires that the ``graphviz`` program is installed, so that the ``dot`` program is in the path. For example, on Debian systems, this can be achieved by typing ``sudo apt install graphviz``, and verified by typing ``dot -V``.

.. note::

    Implementing custom lattices with more than one site per unit cell requires an advanced install.

Testing
~~~~~~~

You can confirm that HofstadterTools is correctly installed by running the pytests.

.. code:: console

    user@domain:any/path$ pytest --pyargs HT

Once the *project* ``HofstadterTools`` is installed, the *package* ``HT`` will be available in your Python environment. In addition, you can access the programs ``band_structure``, ``butterfly``, ``plot_band_structure``, and ``plot_butterfly``, from any location.

.. code:: console

    user@domain:any/path$ band_structure --help
    user@domain:any/path$ butterfly --help
    user@domain:any/path$ plot_band_structure --help
    user@domain:any/path$ plot_butterfly --help

The ``plot_*`` programs are used to replot band_structures / butterflies that have been saved to file.

Directory Structure
-------------------

* **src** -- sources root with the ``HT`` package, along with its configuration settings, subpackages, and programs. A detailed description of the available `programs <https://hofstadter.tools/tutorials.html>`__ and `namespace packages <https://hofstadter.tools/_autosummary/functions.html>`__ is in the documentation.

  * **HT** -- ``HT`` package.

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

Contributions are always welcome! The HofstadterTools repository is maintained using `GitHub <https://github.com/HofstadterTools/HofstadterTools>`__. If you would like to contribute, please submit a `pull request <https://github.com/HofstadterTools/HofstadterTools/pulls>`__; if you would like to report an issue or problem, please open an `issue <https://github.com/HofstadterTools/HofstadterTools/issues>`__; and if you need to seek support, please start a `discussion <https://github.com/HofstadterTools/HofstadterTools/discussions>`__. For all other enquires, please contact `Bart Andrews <https://bartandrews.me>`__.
