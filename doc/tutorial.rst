Tutorial
========

This tutorial will run through some very basic examples of how to do a few
interesting things with the :py:mod:`hyperspy_tools.plotting` and
:py:mod:`hyperspy_tools.shifting_lines` packages.

The code in this section is also available as a Jupyter
:download:`notebook <HyperSpy_Tools_Demos.ipynb>`, for interactive evaluation.

Getting Started
+++++++++++++++

#.  Download the tutorial :download:`files <tutorial_files.zip>` and extract
    them to a folder.

#.  Install the dependencies and the :mod:`~hyperspy_tools` module
    following the :doc:`install` instructions.

#.  Fire up a Python (or Jupyter or Notebook) instance, and import the module
    packages (we also set up a bit of formatting code):

    ..  code-block:: python

        >>> import hyperspy.api as hs
        >>> import hyperspy_tools.plotting as htp
        >>> import hyperspy_tools.shifting_lines as hts

        >>> # Set some display preferences:
        >>> import seaborn as sns
        >>> sns.set_style('white',
        ...               rc={"image.cmap": 'cubehelix',
        ...                   'legend.frameon': False,
        ...                   "lines.linewidth": 1})
        >>> sns.set_context('poster', font_scale=1.2)

        >>> import matplotlib.pyplot as plt

Plotting survey images from DigitalMicrograph
+++++++++++++++++++++++++++++++++++++++++++++

One of the more useful methods in the :py:mod:`hyperspy_tools.plotting`
package is :py:meth:`~hyperspy_tools.plotting.plot_dm3_survey_with_markers`.
HyperSpy currently does not load the annotations present on survey images
that are saved by `DigitalMicrograph <http://www.gatan.com/products/tem-
analysis/gatan-microscopy-suite-software>`_. These annotations show where
the beam was located, from where a spectrum image was acquired, and where
the spatial drift correction was performed.
:py:meth:`~hyperspy_tools.plotting.plot_dm3_survey_with_markers` will load
this information from the ``.dm3`` tags, and plot it in a manner consistent
with HyperSpy:

    .. code-block:: python

        >>> htp.plot_dm3_survey_with_markers('survey_image.dm3')

    .. figure:: figures/survey_image.png
       :width: 500 px
       :alt: DigitalMicrograph Survey Image
       :align: center

See the detailed documentation of :py:meth:`~hyperspy_tools.plotting.plot_dm3_survey_with_markers`
for more information about the different options for plotting.

Customizing plot outlines and colorbars
+++++++++++++++++++++++++++++++++++++++

A few methods are included that are useful for plotting images:
:py:meth:`~hyperspy_tools.plotting.add_colored_outlines` will outline images
with colored borders (useful for matching up with plots of spectra), and
:py:meth:`~hyperspy_tools.plotting.add_custom_colorbars` will easily add colorbars
to image plots with ticks at specified locations. A simple demonstration:

    ..  code-block:: python

        >>> eels_sig = hs.load('EELS_signal.hdf5')
            # This spectrum already had a decomposition performed...
        >>> loadings = eels_sig.get_decomposition_loadings()
        >>> hs.plot.plot_images(loadings,
        ...                     axes_decor=None,
        ...                     per_row=3,
        ...                     label=['Loading {}'.format(i) for i in range(3)],
        ...                     colorbar=None)

    .. figure:: figures/loadings_wo_outlines.png
       :width: 100%
       :alt: Plain plot of decomposition loadings
       :align: center

Add the colored outlines:

    .. code-block:: python

        >>> htp.add_colored_outlines(fig=plt.gcf(),
        ...                          signal=eels_sig,
        ...                          num_images=3,
        ...                          border=0,
        ...                          lw=15)

    ..  figure:: figures/loadings_with_outlines.png
        :width: 100%
        :alt: Decomposition loadings with outlines
        :align: center

Add the custom colorbars:

    .. code-block:: python

        >>> # Little helper function to calculate middle of tick_list easily
        >>> def avg_list(i, j):
        >>>     return [i, (i + j)/2, j]

        >>> # add the colorbars
        >>> htp.add_custom_colorbars(fig=plt.gcf(),
        ...                          tick_list=[avg_list(16, 28),
        ...                                     avg_list(-21, 0),
        ...                                     avg_list(0, 12)])

    ..  figure:: figures/loadings_with_colorbars.png
        :width: 100%
        :alt: Decomposition loadings with outlines
        :align: center


Correcting spatial drift
++++++++++++++++++++++++

When collecting a spectrum image across a planar interface, spatial drift
during the acquisition can cause the interface to appear slanted or jagged.
Spatial drift correction methods in the acquisition software help, but are
not always perfect.

The :py:mod:`hyperspy_tools.shifting_lines` package provides a simple
means to correct this drift, using the STEM signal that is collected at the
same time as the spectrum image. The small example presented here will
demonstrate the process:

#.  Load the STEM signal and the spectrum image:

    .. code-block:: python

        >>> stem = hs.load('STEM_signal.dm3')
        >>> eels = hs.load('EELS_signal.hdf5')

#.  The :py:meth:`~hyperspy_tools.shifting_lines.get_shifts_from_area_stem`
    method will extract individual line profiles from the STEM image, find
    the midpoint of the intensity, and report what spatial shift is necessary
    to bring that midpoint to an average coincident point with all the other
    profiles. The method returns ``stem_linescans`` (a list of the extracted
    line profiles) and ``shifts`` (the array of shift values):

    .. code-block:: python

        >>> stem_linescans, shifts = hts.get_shifts_from_area_stem(stem,
        ...                                                        debug=False)

    a.  If the ``debug=True`` option is provided, a scatter plot of the
        measured midpoints will be shown:

        ..  figure:: figures/midpoints_plot.png
            :width: 100%
            :alt: Plot of line scan midpoints
            :align: center

    b.  The ``shifts`` array now contains all of the shifts necessary to bring
        all the scans together:

        ..  code-block:: python

            >>> shifts
                array([-0.2715, -0.4605, -0.1705, -0.3105, -0.3635, -0.0135, -0.1235,
                       -0.3245, -0.3675,  0.0135,  0.4585,  0.3345,  0.1745, -0.5165,
                        0.0515, -0.1095, -0.2035,  0.3415,  0.5375,  0.2265,  0.3655,
                        0.5255, -0.2275,  0.4275,  0.7455,  0.6245,  0.4805, -0.1155,
                       -0.3335,  0.6915])

#.  The STEM signal can be easily shifted to become planar with the
    :py:meth:`~hyperspy_tools.shifting_lines.shift_area_stem` method.
    HyperSpy's :py:func:`~hyperspy.drawing.utils.plot_images` method can
    be used to visualize the results:

        ..  code-block:: python

            >>> shifted_stem = hts.shift_area_stem(stem, shifts=shifts)

            >>> hs.plot.plot_images([stem, shifted_stem],
            ...                     suptitle='Shifting STEM images',
            ...                     colorbar=None,
            ...                     label=['Original', 'Shifted'],
            ...                     axes_decor=None,
            ...                     scalebar='all')

        ..  figure:: figures/shifted_stem.png
            :width: 100%
            :alt: Comparison of shifted STEM image to original
            :align: center

    As can be seen, the interface is straightened, and the image was cropped
    to prevent the presence of empty pixels in the image. If desired, the
    ``crop_scan`` parameter can be set to ``False``, and the image will not be
    cropped.

#.  An analogous method :py:meth:`~hyperspy_tools.shifting_lines.shift_area_eels`
    allows for a similar operation on EELS spectrum images as well. The
    HyperSpy :py:func:`~hyperspy.drawing.utils.plot_signals` function can
    compare the results:

        ..  code-block:: python

            >>> # Do not crop the scan so we can see the empty pixels:
            >>> shifted_eels = hts.shift_area_eels(eels,
            ...                                    shifts=shifts,
            ...                                    crop_scan=False)

            >>> # Plotting as images using the as_image() transformation
            >>> # (rather than Spectra) makes it easier to see the shift
            >>> hs.plot.plot_signals([eels.as_image((0,1)),
            ...                       shifted_eels.as_image((0,1))])

        ..  figure:: figures/shifted_eels.png
            :width: 100%
            :alt: Comparison of shifted EELS signals to original
            :align: center