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

import hyperspy.api as _hs
import numpy as _np
from tqdm import tqdm as _tqdm
import math as _math
import matplotlib.pyplot as _plt

__all__ = ['extract_ZLP',
           'align_energy_vertical']


def extract_ZLP(signal,
                range_factor=(-10, 6),
                return_inelastic=True,
                print_output=True):
    """
    Extracts the ZLP of an EELS spectrum and returns it and the inelastic
    spectrum.

    ZLP is fit with the sum of a Gaussian and Lorentzian function. Fit
    range is set by first fitting a gaussian and calculating HWHM based
    off that. This is the method used by Gatan's DigitalMicrograph software


    Parameters
    ----------
    signal: ~hyperspy.signal.Signal
        Low-loss signal from which to extract the zero-loss

    range_factor: tuple
        This parameter must be a tuple of length two, with the first value
        less than the second. The values here correspond to the range
        over which the ZLP will be fit. The HWHM of the initial ZLP fit
        will be multiplied by these factors to obtain the fit range.

    return_inelastic: bool
        This switch controls whether the method will return the
        inelastic signal as well, which is calculated by subtracting the
        ZLP from the original signal and then "cleaning up" by removing
        negative values

    print_output: bool
        This switch controls whether or not the method will output any
        results, besides the return.

    Returns
    -------
    zlp: Signal instance representing the ZLP of the spectrum
    inelastic: Signal instant representing the inelastic portion of
        the spectrum (only returned if ``return_inelastic`` is set)

    """
    # check range factor
    if len(range_factor) is not 2 or range_factor[0] >= range_factor[1]:
        print("Range factors were not valid, resetting to default "
              "values of (-10, +6)...")
        range_factor = (-10, 6)

    orig_progress = _hs.preferences.General.show_progressbar

    if print_output:
        _hs.preferences.General.show_progressbar = True
    else:
        _hs.preferences.General.show_progressbar = False

    # lambda function to calculate HWHM from gaussian sigma parameter
    def sigtohwhm(sig):
        return 2 * _np.sqrt(2 * _np.log(2)) * sig / 2
    m1 = signal.create_model()  # create model of signal
    G = _hs.model.components1D.Gaussian()  # create gaussian component
    # add gaussian component to model
    m1.append(G)
    # deactivate the power law background
    m1['PowerLaw'].active = False

    # set fitting range to +-6 eV for initial fitting
    m1.set_signal_range(-6, 6)
    if print_output:
        print("Fitting initial Gaussian to determine fitting range")
        # perform initial gaussian fit
        m1.multifit(show_progressbar=True)
    else:
        m1.multifit(show_progressbar=False)

    # copy fit gaussian's sigma into new signal
    HWHM = G.sigma.as_signal()
    # save the centre value of G for later (helps pre-fit total ZLP)
    G_centre = G.centre.as_signal()
    # convert to HWHM from sigma by mapping
    HWHM.map(sigtohwhm, show_progressbar=False)
    HWHM.metadata.General.title = 'HWHM of Gaussian component'
    # calculate mean of HWHM across signal
    HWHM_mean = HWHM.data.mean()
    if print_output:
        print("Mean HWHM of ZLP is {0} {1}".format(
            round(HWHM_mean, 2), signal.axes_manager[-1].units))

    # perform actual ZLP fitting
    # create model and add components
    m2 = signal.create_model()
    g1 = _hs.model.components1D.Gaussian()
    l1 = _hs.model.components1D.Lorentzian()
    m2.extend((g1, l1))
    m2['PowerLaw'].active = False
    g1.centre.map['values'][:] = G_centre.data
    g1.centre.map['is_set'][:] = True

    # set signal range (fixed from range_factor)
    m2.set_signal_range(
        range_factor[0] *
        HWHM_mean,
        range_factor[1] *
        HWHM_mean)
    if print_output:
        print("Fitting Signal in range " + repr(range_factor) + "*HWHM...")
        m2.multifit(show_progressbar=True)
        print("Fit results:")
        m2.print_current_values()
        print("")
    else:
        m2.multifit(show_progressbar=False)

    # calculate inelastic part by subtracting modeled ZLP:
    if print_output:
        print("Calculating inelastic portion of spectrum:")

    inelastic = signal - m2.as_signal()
    # set all values in fit range for inelastic to 0; remove any
    # remaining negative residuals
    inelastic.isig[:range_factor[1] * HWHM_mean] = 0
    inelastic.data[inelastic.data < 0] = 0
    inelastic.metadata.General.title += ' - inelastic'  # change tiles

    # calculate ZLP by removing inelastic from original signal
    if print_output:
        print("Calculating ZLP portion of spectrum:")
    zlp = signal - inelastic
    zlp.metadata.General.title += ' - ZLP'

    # calculating "residuals"
    if print_output:
        model = m2.as_signal().isig[range_factor[0] *
                                    HWHM_mean:range_factor[1] * HWHM_mean].data
        real = signal.isig[range_factor[0] *
                           HWHM_mean:range_factor[1] * HWHM_mean].data

        print("\"Residual\" mean: {0}".format(
            round(abs(real - model).mean(), 2)))
        print("\"Residual\" std dev: {0}".format(
            round(abs(real - model).std(), 2)))

    _hs.preferences.General.show_progressbar = orig_progress

    if return_inelastic:
        return zlp, inelastic
    else:
        return zlp


def align_energy_vertical(signal,
                          start=None,
                          end=None,
                          column=0,
                          smoothing_parameter=0.05,
                          number_of_iterations=1,
                          print_output=True,
                          plot_deriv=False,
                          plot_shifts=False):
    """
    Align the energy (signal) axis of a spectrum image, based off of the
    spectrum in the first column of each row. This is useful for SI of
    vertically aligned interfaces.

    If the spectrum image is uniform in the left-most column, this method
    should help to reduce the effect of energy drift that occurs during
    acquisition (if you cannot acquire the zero-loss at the same time). It will
    not correct for any drift that occurs within each row, but that is a
    much trickier problem, considering the SI likely covers different phases.

    A note on the method: the first column is extracted as a line scan. This
    data is smoothed using a Lowess filter (to reduce the impact of noise), and
    then the signal-dimension derivative is taken. This signal is passed to
    hyperspy.signal.Signal1DTools.estimate_shift1D in order to figure out
    what shift is necessary to keep the energy axis aligned. Each column
    of the original spectrum image is shifted by this same amount, and then
    the overall spectrum image is cropped in the signal dimension so there
    are no blank pixels. If the results are not quite as expected,
    try increasing the smoothing parameter, as noise in the derivative is
    the most likely reason for failure.

    Parameters
    ----------
    signal: ~hyperspy.signal.Signal
        2D spectrum image to align in energy dimension
    start: {int | float | None}
        The limits of the interval in which to align. If int they are taken
        as the axis index. If float they are taken as the axis value.
    end: {int | float | None}
        The limits of the interval in which to align. If int they are taken
        as the axis index. If float they are taken as the axis value.
    column: int
        The column of data to use for shifting. By default, the left-most (
        0) is used, but if there is no edge in this area, a different
        column should be used.
    smoothing_parameter: float
        Degree of smoothing used to smooth the original spectral data
        (necessary before taking the derivative). This parameter is passed to
        :py:meth:`~hyperspy.signal.Signal1DTools.smooth_lowess`
    number_of_iterations: int
        Number of Lowess iterations used to smooth the original spectral data
        (necessary before taking the derivative). This parameter is passed to
        :py:meth:`~hyperspy.signal.Signal1DTools.smooth_lowess`
    print_output: bool
        Whether or not to show output during calculation.
    plot_deriv: bool
        Whether or not to plot the derivative output. Useful if results are
        not as expected, and can show if more smoothing is needed
    plot_shifts: bool
        Whether or not to show a plot illustrating the shifts that were found

    Returns
    -------
    aligned_signal: ~hyperspy.signal.Signal
        2D spectrum image with signal axes aligned and cropped
    """
    s = signal.inav[column, :]
    if print_output:
        print("Smoothing column {}...".format(column))
    s.smooth_lowess(smoothing_parameter=smoothing_parameter,
                    number_of_iterations=number_of_iterations,
                    show_progressbar=True)
    sd = s.diff(-1)
    if plot_deriv:
        sd.plot()

    shifts = sd.estimate_shift1D(start=start, end=end)

    if print_output:
        print('Shifts are:')
        print(shifts)
        print('Max shift: {}'.format(_np.nanmax(shifts)))

    if plot_deriv:
        sdd = sd.deepcopy()
        sdd.shift1D(shifts, crop=False, show_progressbar=False)
        sdd.plot()

    if plot_shifts:
        u = s.axes_manager[-1].units
        med = _np.median(shifts)
        _plt.figure()
        _plt.scatter(range(len(shifts)), shifts)
        _plt.axhline(med, ls='--', c='k')
        _plt.text(0.5, med + 0.05, 'Median = {0:.2f} {1:s}'.format(med, u))
        ax = _plt.gca()
        ax.set_ylabel('Shift value({:s})'.format(u))
        ax.set_xlabel('Row #')
        _plt.xlim(0, len(shifts))
        print(("Median shift is {:.2f}".format(_np.median(shifts))))
        print(("Mean shift is {:.2f}".format(_np.mean(shifts))))

    aligned_signal = signal.deepcopy()

    for i in _tqdm(range(signal.axes_manager['x'].size), desc='Aligning '
                                                              'spectrum '
                                                              'image'):
        s = signal.inav[i, :]
        s.shift1D(shifts, crop=False, show_progressbar=False)
        aligned_signal.inav[i, :] = s

    # ## This code lifted from HyperSpy's shift1D method:
    # Figure out min/max shifts, and translate to shifts in index as well
    minimum, maximum = _np.nanmin(shifts), _np.nanmax(shifts)
    axis = aligned_signal.axes_manager.signal_axes[0]
    if minimum < 0:
        ihigh = 1 + axis.value2index(
            axis.high_value + minimum,
            rounding=_math.floor)
    else:
        ihigh = axis.high_index + 1
    if maximum > 0:
        ilow = axis.value2index(axis.offset + maximum,
                                rounding=_math.ceil)
    else:
        ilow = axis.low_index

    aligned_signal.crop(axis.index_in_axes_manager,
                        ilow,
                        ihigh)

    return aligned_signal

