#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# HyperSpy Tools documentation build configuration file, created by
# sphinx-quickstart on Thu Apr  7 11:40:21 2016.

import sys
import os
import io
from os import path
import re
import sphinx_rtd_theme
from mock import Mock

MOCK_MODULES = ['numpy', 'sip', 'matplotlib', 'matplotlib.pyplot',
                'mpl_toolkits', 'mpl_toolkits.axes_grid1',
                'matplotlib.patches', 'hyperspy', 'hyperspy.api',
                'hyperspy.io_plugins',
                'hyperspy.io_plugins.digital_micrograph', 'PyQt4',
                'PyQt4.QtCore', 'PyQt4.QtGui', 'seaborn', 'hyperspy', 'scipy',
                'scipy.interpolate']

for mod_name in MOCK_MODULES:
    sys.modules[mod_name] = Mock()

sys.path.insert(0, os.path.abspath('..'))


def setup(app):
    app.add_javascript('copybutton.js')


# next two methods read version from file
def read(*names, **kwargs):
    with io.open(
            path.join(path.dirname(__file__), *names),
            encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


# ################ EDIT THIS ############################
version = find_version('../hyperspy_tools/__init__.py')
release = version
# ################ EDIT THIS ############################

intersphinx_mapping = {'scipy': ('http://docs.scipy.org/doc/scipy/reference/',
                                 None),
                       'hyperspy': ('http://hyperspy.org/hyperspy-doc/dev/',
                                    None),
                       'python': ('http://docs.python.org/3.5', None),
                       'numpy': ('http://docs.scipy.org/doc/numpy-1.10.0',
                                 None),
                       'matplotlib': ('http://matplotlib.org/', None),
                       'seaborn': (
                           'https://stanford.edu/~mwaskom/software/seaborn'
                           '/', None)}

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinxcontrib.napoleon',
]

templates_path = ['_templates']
html_static_path = ['_static']
html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
source_suffix = '.rst'


# The master toctree document.
master_doc = 'index'

# General information about the project.
project = 'HyperSpy Tools'
copyright = '2016, Joshua Taillon'
author = 'Joshua Taillon'


language = None
exclude_patterns = ['_build']
pygments_style = 'sphinx'

todo_include_todos = False

htmlhelp_basename = 'HyperSpyToolsdoc'

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
}

latex_documents = [
    (master_doc, 'HyperSpyTools.tex', 'HyperSpy Tools Documentation',
     'Joshua Taillon', 'manual'),
]


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'hyperspytools', 'HyperSpy Tools Documentation',
     [author], 1)
]

# If true, show URL addresses after external links.
#man_show_urls = False


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'HyperSpyTools', 'HyperSpy Tools Documentation',
     author, 'HyperSpyTools', 'One line description of project.',
     'Miscellaneous'),
]
