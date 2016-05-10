Introduction
============
This module is used to perform some various tasks that are useful to me when
using `HyperSpy`_, but probably not general enough for inclusion in the primary
codebase. The package is currently split into three different modules:
:py:mod:`~hyperspy_tools.plotting` provides some useful extra plotting
configuration; :py:mod:`~hyperspy_tools.shifting_lines` contains
code that can be used to align an area scan (or line scans) across a planar
interface that appear to jump around due to spatial drift in the TEM;
:py:mod:`~hyperspy_tools.eels` provides a few methods useful when dealing
with EELS data, beyond what is available in HyperSpy by default, including
a method to correct energy drift of spectrum images and to extract zero-loss
data.


.. _HyperSpy: http://hyperspy.org/hyperspy-doc/dev/
