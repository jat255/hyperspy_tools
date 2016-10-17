Installation
============

Requirements
++++++++++++

The following dependencies are needed to run the code (as-written). Some of
these could be removed, but the code was written with my research
results in mind, and as such, has some dependencies that make things match
my personal preferences (such as ``seaborn``). Some details on installing them
are given below:

..  csv-table::
    :header: Package, Source, PyPI
    :escape: \

    ``numpy``, `Source <Numpy_>`_, `PyPI <NumpyPYPI_>`_
    ``matplotlib``, `Source <matplotlib_>`_, `PyPI <matplotlibPYPI_>`_
    ``scipy``, `Source <scipy_>`_, `PyPI <scipyPYPI_>`_
    ``hyperspy``, `Source <HyperSpy_>`_, `PyPI <HyperSpyPYPI_>`_
    ``seaborn``, `Source <seaborn_>`_, `PyPI <seabornPYPI_>`_
    ``tqdm``, `Source <tqdm_>`_, `PyPI <tqdmPYPI_>`_
    ``PyQt4``, `Source <PyQt4_>`_,
    ``h5py``, `Source <h5py_>`_, `PyPI <h5pyPYPI_>`_
    ``scikit-image``, `Source <skimage_>`_ , `PyPI <skimagePYPI_>`_
    ``docopt``, `Source <docopt_>`_, `PyPI <docoptPYPI_>`_
    ``pandas``, `Source <pandas_>`_, `PyPI <pandasPYPI_>`_

.. _Numpy: http://www.numpy.org/
.. _NumpyPYPI: https://pypi.python.org/pypi/numpy/1.11.0

.. _matplotlib: http://matplotlib.org/
.. _matplotlibPYPI: https://pypi.python.org/pypi/matplotlib/1.5.1

.. _scipy: https://www.scipy.org/scipylib/index.html
.. _scipyPYPI: https://pypi.python.org/pypi/scipy/0.7.0

.. _HyperSpy: http://www.hyperspy.org/
.. _HyperSpyPYPI: https://pypi.python.org/pypi/hyperspy/0.8.4

.. _seaborn: https://stanford.edu/~mwaskom/software/seaborn/
.. _seabornPYPI: https://pypi.python.org/pypi/seaborn

.. _tqdm: https://github.com/tqdm/tqdm/
.. _tqdmPYPI: https://pypi.python.org/pypi/tqdm

.. _PyQt4: https://www.riverbankcomputing.com/software/pyqt/download

.. _h5py: http://www.h5py.org/
.. _h5pyPYPI: https://pypi.python.org/pypi/h5py/2.5.0

.. _skimage: https://github.com/scikit-image/scikit-image
.. _skimagePYPI: https://pypi.python.org/pypi/scikit-image

.. _docopt: https://github.com/docopt/docopt
.. _docoptPYPI: https://pypi.python.org/pypi/docopt/0.6.2

.. _pandas: http://pandas.pydata.org/
.. _pandasPYPI: https://pypi.python.org/pypi/pandas/0.18.0

.. _statsmodels: http://statsmodels.sourceforge.net/
.. _statsmodelsPYPI: https://pypi.python.org/pypi/statsmodels

.. _OpenCV: http://opencv.org/


Development Version Installation
++++++++++++++++++++++++++++++++

Currently, only the development version is capable of installation.
The latest version of the code should be available in the online
`repository <https://github.com/jat255/hyperspy_tools>`_.
To get this version installed on your system: install the requirements from
above, clone the repository, and then install with ``pip``:

.. code-block:: bash

    $ git clone https://github.com/jat255/hyperspy_tools.git
    $ cd hyperspy_tools
    $ pip install -e ./

