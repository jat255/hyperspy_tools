# Copyright 2015 Joshua Taillon
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from setuptools import setup, find_packages  # Always prefer setuptools
from codecs import open  # To use a consistent encoding
from os import path
import io
import re

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()


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

setup(
    name='hyperspy_tools',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=find_version('hyperspy_tools/__init__.py'),

    description='Python modules providing useful tools for use in HyperSpy',
    long_description='This package supplies a reference point for useful '
                     'scripts that have been developed for use with '
                     'HyperSpy, a hyperspectral data analysis package. '
                     '(see http://www.hyperspy.org)',

    # The project's main homepage.
    url='https://github.com/jat255/hyperspy_tools',

    # Author details
    author='Joshua Taillon',
    author_email='jat255@gmail.com',

    # Choose your license
    license='GPLv3',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: Apache Software License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 2.7',
    ],

    # What does your project relate to?
    keywords='electronmicroscopy TEM EELS',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),

    # List run-time dependencies here. These will be installed by pip when your
    # project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=[],

    extras_require={
        'full': ['numpy >= 1.9.1', #
                 'matplotlib>=1.4.2', #
                 'scipy>=0.15.1', #
                 'hyperspy>=0.8', #
                 'seaborn', #
                 'tqdm', #
                 'PyQt4>=4.11.3',
                 'h5py',
                 'scikit-image',
                 'docopt',
                 'pandas'], #
        'base': [],
    },
)
