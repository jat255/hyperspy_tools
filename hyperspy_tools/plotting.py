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
#
# #########################################################################
# This file contains the code necessary to load a DM survey image and plot
# it with Hyperspy, adding markers as necessary to show the spatial extents
# of the spectrum image, beam location, and spatial drift box
# #########################################################################

import seaborn as _sns
import matplotlib
from matplotlib.patches import Rectangle as _Rect
from mpl_toolkits.axes_grid1 import make_axes_locatable as _mal
import hyperspy

import hyperspy.api as _hs
import hyperspy.io_plugins.digital_micrograph as _dm


__all__ = ['fit_peak',
           'add_colored_outlines',
           'add_custom_colorbars',
           'plot_dm3_survey_with_markers']


def fit_peak(sig, lower_bound, upper_bound, factor_num=None,):
    """
    Fit a Gaussian peak in a range to a signal

    Parameters
    ----------
    sig: ~hyperspy.signal.Signal
        Signal to fit
    lower_bound: float
        Lower bound of values in which to fit
    upper_bound: float
        Upper bound of values in which to fit
    factor_num: int
        if given, the fit will be performed on a decomposition component of
        the signal (given by ``factor_num``), rather than the signal itself

    Returns
    -------
    centre: float
        center of the fitted Gaussian
    """
    if factor_num is not None:
        c1 = sig.get_decomposition_factors()[factor_num]
    else:
        c1 = sig

    # noinspection PyUnresolvedReferences,PyProtectedMember
    if isinstance(sig, hyperspy._signals.spectrum.Spectrum):
        is_eels = False
    else:
        is_eels = True

    if is_eels:
        c1.set_microscope_parameters(beam_energy=200,
                                     convergence_angle=12,
                                     collection_angle=29)
        m1 = c1.create_model(auto_background=False)
    else:
        m1 = c1.create_model()
    g1 = _hs.model.components.Gaussian(centre=((float(lower_bound) +
                                                upper_bound) / 2.0))
    m1.append(g1)
    m1.set_signal_range(lower_bound, upper_bound)
    m1.fit()
    m1.plot()

    centre = g1.centre.value

    return centre


def add_colored_outlines(fig,
                         signal,
                         num_images,
                         color_palette=_sns.color_palette(),
                         border=0.0,
                         lw=15):
    """
    Add outlines to a matplotlib figure to make it easy to visualize with
    the plotted spectra

    Parameters
    ----------
    fig: matplotlib.figure
        figure to which to add outlines (this should not have colorbars,
        add them later)
    signal: ~hyperspy.signal.Signal
        signal that has calibrated axes (in order to set bounds of rectangle)
    num_images: int
        number of images in figure
    color_palette: list
        list of colors to use for outlines
    border: float
        offset of rectangle from edge of data
    lw: float
        width of rectangle lines to plot
    """
    x, y = signal.axes_manager[0].high_value, signal.axes_manager[1].high_value

    for i in range(num_images):
        ax = fig.get_axes()[num_images - (1 + i)]
        # noinspection PyUnresolvedReferences
        r = _Rect((border, border),
                  x - 1.5 * border,
                  y - 1.5 * border,
                  fill=False,
                  edgecolor=color_palette[num_images - (1 + i)],
                  alpha=1,
                  linewidth=lw)
        ax.add_patch(r)


def add_custom_colorbars(fig,
                         tick_list=None):
    """
    Add custom colorbars with ticks at specified values

    Parameters
    ----------
    fig: matplotlib.figure
        figure to which to add colorbars (should not currently have colorbars)
    tick_list: list
        nested list with the position of ticks to be added to colorbars;
        should have a length equal to the number of images in the figure
        Example for a four plot figure:
            >>> tick_list =  [[120,200,280],
            ...               [4,20,36],
            ...               [0,8,16],
            ...               [0,22,44]]
    """
    for i, a in enumerate(fig.axes):
        # if i == 2:
        #     a.get_images()[0].set_clim(0,16)
        div = _mal(a)
        cax = div.append_axes("right", size="5%", pad=0.05)
        if tick_list is None:
            _ = fig.colorbar(a.get_images()[0], cax=cax)
        else:
            # noinspection PyUnresolvedReferences
            _ = fig.colorbar(a.get_images()[0], cax=cax, ticks=tick_list[i])


def plot_dm3_survey_with_markers(fname,
                                 add_text=True,
                                 plot_beam=True,
                                 plot_si=True,
                                 plot_drift=True,
                                 x_offset=1.0,
                                 y_offset=1.0,
                                 im_scale=1.0,
                                 text_size='xx-small',
                                 **kwargs):
    """
    Plot a hyperspy signal with the markers from digital micrograph enabled

    Parameters
    ----------
    fname : str
        Name of .dm3 file to load
    add_text : bool
        Switch to control if labels for the markers are added to the plot
    plot_beam : bool
        Switch to control if beam point is plotted (if present)
    plot_si : bool
        Switch to control if spectrum image box or line is plotted (if present)
    plot_drift : bool
        Switch to control if spatial drift box is plotted (if present)
    x_offset : float
        multiplier to control how far the text will be offset from its
        default position (in the x-direction)
    y_offset : float
        multiplier to control how far the text will be offset from its
        default position (in the y-direction)
    im_scale : float
        will scale the survey image by a given factor (useful for correcting
        image scale bars if the calibration is incorrect)
    text_size : str or float
        size of the text that will be written on the image (follows same
        convention as the `Text
        <http://matplotlib.org/1.3.0/api/artist_api.html#matplotlib.text.
        Text.set_size>`_ matplotlib Artist
    **kwargs
        Other keyword arguments are passed to hs.signal.plot()

    Returns
    -------
    im: hyperspy._signals.image.Image
        HyperSpy signal with markers added. While the figure is open,
        the figure can be saved with a call such as
        ``im._plot.signal_plot.figure.savefig()``
    """

    def _add_beam(image, location, annotationtype):
        if plot_beam:
            beam_m = _hs.plot.markers.point(x=location[1],
                                            y=location[0],
                                            color='red')
            image.add_marker(beam_m)
            if add_text:
                beam_text_m = _hs.plot.markers.text(x=location[1] - (0.5 *
                                                                     x_offset),
                                                    y=location[0] - (1.5 *
                                                                    y_offset),
                                                    color='red',
                                                    text='Beam',
                                                    size=text_size)
                image.add_marker(beam_text_m)
        else:
            pass

    def _add_si(image, location, annotationtype):
        # adds a green rectangle (or line, if the coordinates are such) to
        # image
        if plot_si:
            if annotationtype == 23: # Map
                si_m = _hs.plot.markers.rectangle(x1=location[1],
                                                 y1=location[0],
                                                 x2=location[3],
                                                 y2=location[2],
                                                 color='#13FF00')
            elif annotationtype == 25: # Line profile
                si_m = _hs.plot.markers.line_segment(x1=location[1],
                                                    y1=location[0],
                                                    x2=location[3],
                                                    y2=location[2],
                                                    color='#13FF00')
            else:
                raise Exception("No spectrum image annotation found")
            image.add_marker(si_m)
            if add_text:
                si_text_m = _hs.plot.markers.text(x=location[1],
                                                 y=location[0] - (0.5 *
                                                                  y_offset),
                                                 color='#13FF00',
                                                 text='Spectrum Image',
                                                 size=text_size)
                image.add_marker(si_text_m)
        else:
            pass

    def _add_drift(image, location, annotationtype):
        if plot_drift:
            drift_m = _hs.plot.markers.rectangle(x1=location[1],
                                                 y1=location[0],
                                                 x2=location[3],
                                                 y2=location[2],
                                                 color='yellow')
            image.add_marker(drift_m)
            if add_text:
                drift_text_m = _hs.plot.markers.text(x=location[1],
                                                     y=location[0] - (0.5 *
                                                                     y_offset),
                                                     color='yellow',
                                                     text='Spatial Drift',
                                                     size=text_size)
                image.add_marker(drift_text_m)
        else:
            pass

    im = _hs.load(fname)
    flist = _dm.file_reader(fname)
    annotation_list = flist[0]['original_metadata']['DocumentObjectList'][
                               'TagGroup0']['AnnotationGroupList']

    im.axes_manager[0].scale *= im_scale
    im.axes_manager[1].scale *= im_scale
    scale = im.axes_manager[0].scale

    mapping = {
        'Beam': _add_beam,
        'Spectrum Image': _add_si,
        'Spatial Drift': _add_drift
    }

    im.plot(cmap='gist_gray', **kwargs)

    for i in range(len(annotation_list)):
        try:
            label = annotation_list['TagGroup' + str(i)]['Label']
            loc = annotation_list['TagGroup' + str(i)]['Rectangle']
            scaled_loc = [scale * i for i in loc]
            annotationtype = annotation_list['TagGroup' + str(i)]['AnnotationType']
            mapping[label](im, scaled_loc, annotationtype)

        except KeyError:
            pass

    return im   # Returns a hyperspy._signals.image.Image, and with it all
    #           editing and saving capabilities.


def plot_decomposition_results(signal,
                               list_of_comps,
                               labels=None,
                               factor_style='cascade',
                               loadings_cmap='cubehelix',
                               mva_type='decomposition',
                               loadings_kwargs={},
                               factors_kwargs={}):
    """
    Plot the results of a area SI decomposition (components and factors)

    Parameters
    ----------
    signal: ~hyperspy.signal.Signal
        2D signal for which to plot decomposition results
    list_of_comps: list
        should be a list of ints that defines which factors to plot. Can also
        have a tuple (or another iterable) in place of an int, in which case
        the prescribed components will be summed in the output
    labels: list
        list of labels to use for components
    factor_style: string
        string defining which plot_spectra style to use
    loadings_cmap: string
        colormap to use for the loadings plots
    mva_type: string
        'decomposition' or 'bss'. Defines which HyperSpy method to use
    loadings_kwargs: dict
        other keyword arguments passed to plot_images for the loadings plot
    factors_kwargs: dict
        other keyword arguments passed to plot_spectra for the factors plot
    """

    if mva_type == 'decomposition':
        l = signal.get_decomposition_loadings()
        f = signal.get_decomposition_factors()
    elif mva_type == 'bss':
        l = signal.get_bss_loadings()
        f = signal.get_bss_factors()
    else:
        raise ValueError('Did not understand mva_type: {}'.format(mva_type))

    l_to_plot = []
    f_to_plot = []
    for i in list_of_comps:
        try:
            l_x = [l.inav[j] for j in iter(i)]
            f_x = [f.inav[j] for j in iter(i)]
            summed_comp = l_x[0]
            summed_fact = f_x[0]
            for j in range(1, len(l_x)):
                summed_comp += l_x[j]
                summed_fact += f_x[j]
            l_to_plot.append(summed_comp)
            f_to_plot.append(summed_fact)
        except:
            l_to_plot.append(l.inav[i])
            f_to_plot.append(f.inav[i])

    _hs.plot.plot_images(l_to_plot,
                         cmap=loadings_cmap,
                         label=labels,
                         **loadings_kwargs)

    _hs.plot.plot_spectra(f_to_plot,
                          style=factor_style,
                          legend=labels,
                          **factors_kwargs)
