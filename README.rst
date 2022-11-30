HofstadterTools
===============

HofstadterTools is a numerical toolbox for analyzing the Hofstadter model.

Python environment
------------------

In order to use the tried-and-tested packages for HofstadterTools, we recommend the conda environment ``HT`` in ``environment.yml``.

1) If you have not already, you need to install Anaconda or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`__ (recommended).

2) On the first use, you need to create the ``HT`` environment: ``conda env create -f environment.yml``

3) Whenever you would like to use the environment, you can then run: ``conda activate HT``

We only include packages in the environment that are necessary for HofstadterTools. Currently installed packages:

- python
- numpy
- matplotlib
- prettytable
- tqdm

- furo
- ipython
- sphinx
- sphinx_rtd_theme
- pydata-sphinx-theme
- sphinx-autodoc-typehints
- nbsphinx
- jupyter
- jupytext

References
----------

`[Bauer2022] <https://arxiv.org/abs/2110.09565>`__ "Fractional Chern insulators with a non-Landau level continuum limit", by David Bauer et al., PRB **105**, 045144 (2022).
