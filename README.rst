HofstadterTools
===============

**H**\ ofstadter\ **T**\ ools (\ **HT**) is a set of python programs and classes for analyzing the Hofstadter model.

Quick Start
-----------

Using HofstadterTools is quick and easy! Assuming the recommended scenario of a UNIX shell with `git <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`__ and `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`__ installed:

.. code:: console

		(base) user@domain:~$ git clone https://github.com/bartandrews/HofstadterTools.git
		(base) user@domain:~$ cd HofstadterTools
		(base) user@domain:~/HofstadterTools$ conda env create -f environment.yml
		(base) user@domain:~/HofstadterTools$ conda activate HT
		(HT) user@domain:~/HofstadterTools$ cd code
		(HT) user@domain:~/HofstadterTools/code$ python band_structure.py -mod Hofstadter -nphi 1 4

Voilà! You have just plotted the band structure of the Hofstadter model with a flux density of :math:`n_\phi=1/4`. Now you can explore the :doc:`programs <programs>` and :ref:`code reference <code_reference>` to see what HofstadterTools has to offer.

Python Environment
------------------

We recommend the use of a python virtual environment to handle the package dependencies. In the following, we assume a UNIX shell, however these instructions may be readily adapted for Windows.

Using ``conda``:

1) If you have not already, install Anaconda or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`__ (recommended).
2) On first use, create the ``HT`` environment: ``conda env create -f environment.yml``
3) Whenever you would like to use the environment, run: ``conda activate HT``

Using ``pip``:

1) On most UNIX-based operating systems, ``pip`` is already installed. If not, install `pip <https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#installing-pip>`__.
2) Create the virtual environment (recommended in the project root): ``python -m venv env``
3) Activate the virtual environment: ``source env/bin/activate``
4) Install the dependencies: ``pip install -r requirements.txt``

.. note::

		For compiling the documentation, ``sphinx_rtd_theme`` was installed using pip to get the newer version number (>=0.5.1). This fixes a minor bug with the formatting of unordered lists.
