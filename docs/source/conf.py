# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from datetime import datetime
# import inspect
sys.path.insert(0, os.path.abspath('../../src/HT'))  # Source code dir relative to this file

# -- Project information -----------------------------------------------------------------------------------------------

project = 'HofstadterTools'
copyright = f'{datetime.now().year}, Bartholomew Andrews'
author = 'Bartholomew Andrews'

# The full version, including alpha/beta/rc tags
release = '1.0.0'

# -- pandoc ------------------------------------------------------------------------------------------------------------

from inspect import getsourcefile

# Get path to directory containing this file, conf.py.
DOCS_DIRECTORY = os.path.dirname(os.path.abspath(getsourcefile(lambda: 0)))


def ensure_pandoc_installed(_):
    import pypandoc

    # Download pandoc if necessary. If pandoc is already installed and on
    # the PATH, the installed version will be used. Otherwise, we will
    # download a copy of pandoc into docs/bin/ and add that to our PATH.
    pandoc_dir = os.path.join(DOCS_DIRECTORY, "bin")
    # Add dir containing pandoc binary to the PATH environment variable
    if pandoc_dir not in os.environ["PATH"].split(os.pathsep):
        os.environ["PATH"] += os.pathsep + pandoc_dir
    pypandoc.ensure_pandoc_installed(
        targetfolder=pandoc_dir,
        delete_installer=True,
    )


def setup(app):
    app.connect("builder-inited", ensure_pandoc_installed)

# -- General configuration ---------------------------------------------------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',  # Core library for html generation from docstrings
    'sphinx.ext.autosummary',  # Create neat summary tables
    'sphinx.ext.doctest',
    'sphinx.ext.extlinks',
    'sphinx.ext.intersphinx',  # Link to other project's documentation (see mapping below)
    'sphinx.ext.todo',
    'sphinx.ext.mathjax',
    'sphinx.ext.graphviz',
    'sphinx.ext.inheritance_diagram',
    'sphinx_rtd_theme',
    'sphinx_copybutton',
    'sphinxcontrib.bibtex',
    'sphinx.ext.viewcode',  # Add a link to the Python source code for classes, functions etc.
    'sphinx_autodoc_typehints',  # Automatically document param types (less noise in class signature)
    'nbsphinx',  # Integrate Jupyter Notebooks and Sphinx
    'IPython.sphinxext.ipython_console_highlighting',
    'sphinx.ext.napoleon'
]

autosummary_generate = True  # Turn on sphinx.ext.autosummary
# autoclass_content = "both"  # Add __init__ doc (i.e. params) to class summaries
html_show_sourcelink = False  # Remove 'view source code' from top of page (for html, not python)
autodoc_inherit_docstrings = True  # If no docstring, inherit from base class
set_type_checking_flag = True  # Enable 'expensive' imports for sphinx_autodoc_typehints
nbsphinx_allow_errors = True  # Continue through Jupyter errors
# autodoc_typehints = "description" # Sphinx-native method. Not as good as sphinx_autodoc_typehints
add_module_names = False  # Remove namespaces from class/method signatures

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes. e.g. furo
#
# Readthedocs theme
# on_rtd is whether on readthedocs.org, this line of code grabbed from docs.readthedocs.org...
on_rtd = os.environ.get("READTHEDOCS", None) == "True"
if not on_rtd:  # only import and set the theme if we're building docs locally
    import sphinx_rtd_theme
    html_theme = "sphinx_rtd_theme"
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
    # html_theme = "furo"

html_logo = "images/logo.png"
html_favicon = "images/logo.png"

html_css_files = ["custom.css"]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_context = {
    "display_github": True,  # Integrate GitHub
    "github_user": "bartandrews",  # Username
    "github_repo": "HofstadterTools",  # Repo name
    "github_version": "main",  # Version
    "conf_py_path": "/docs/source/",  # Path in the checkout to the docs root
}

html_theme_options = {
    'collapse_navigation': False,
    'style_external_links': True,
}

# EPUB options
epub_show_urls = 'footnote'

# -- sphinx.ext.intersphinx --------------------------------------------------------------------------------------------
# cross links to other sphinx documentations
# this makes  e.g. :class:`numpy.ndarray` work

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference/', None),
    'matplotlib': ('https://matplotlib.org', None),
    'h5py': ('https://docs.h5py.org/en/stable/', None),
}

# -- sphinxcontrib.bibtex ----------------------------------------------------------------------------------------------

bibtex_bibfiles = ['references.bib']

from pybtex.style.formatting.unsrt import Style as UnsrtStyle
from pybtex.style.labels import BaseLabelStyle
from pybtex.style.sorting.author_year_title import SortingStyle
from pybtex.plugin import register_plugin


class CustomBibtexStyle1(UnsrtStyle):
    default_sorting_style = 'key'
    default_label_style = 'key'


class CustomBibtexStyle2(UnsrtStyle):
    default_sorting_style = 'year_author_title'
    default_label_style = 'key'


class KeyLabelStyle(BaseLabelStyle):
    def format_labels(self, sorted_entries):
        return [entry.key for entry in sorted_entries]


class YearAuthorTitleSort(SortingStyle):
    def sorting_key(self, entry):
        author_key, year, title = super().sorting_key(entry)
        return year, author_key, title


class KeySort(SortingStyle):
    def sorting_key(self, entry):
        return entry.key


register_plugin('pybtex.style.formatting', 'custom1', CustomBibtexStyle1)
register_plugin('pybtex.style.formatting', 'custom2', CustomBibtexStyle2)
register_plugin('pybtex.style.labels', 'key', KeyLabelStyle)
register_plugin('pybtex.style.sorting', 'key', KeySort)
register_plugin('pybtex.style.sorting', 'year_author_title', YearAuthorTitleSort)

# -- sphinx.ext.inheritance_diagram ------------------------------------------------------------------------------------

inheritance_graph_attrs = {
    'rankdir': "TB",  # top-to-bottom
    'fontsize': 14,
    'ratio': 'compress',
}

# -- sphinx.ext.napoleon -----------------------------------------------------------------------------------------------
# numpy-like doc strings

napoleon_use_admonition_for_examples = True
napoleon_use_ivar = False  # otherwise :attr:`...` doesn't work anymore
napoleon_custom_sections = ['Options']
napoleon_attr_attributes = True
