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
# ##################################################################
# This file contains the code necessary to measure the spatial shift
# of EELS line scans from their corresponding STEM signal, and also
# to shift the line scans and then build them into an array like an
# areal spectrum image.
# #################################################################

import numpy as np
import sys
import re
from pprint import pprint
from scipy import stats
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import hyperspy.api as hs
from PyQt4 import QtGui, QtCore
from progressbar import ProgressBar, Percentage, Bar, ETA


def _natural_sort(l):
    """
    Internal function to sort a list naturally (10 comes after 9)
    (from: http://stackoverflow.com/questions/4836710/does-python-have-a-
    built-in-function-for-string-natural-sort)

    Parameters:
    -----------
    l: list
        list of values to sort naturally

    Returns:
    --------
    sorted list
    """
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


# noinspection PyTypeCheckerInspection
def gui_fnames(title=None, directory=None):
    """
    Select a file via a dialog and returns the file name.
    If multiple selected, returns list of file names.

    Parameters
    ----------
    title: str or None
        String to use as title of file chooser window
    directory: str or None
        Directory to use as the starting directory for the file chooser

    Returns:
    --------
    fname: str or list
        single string or list of strings representing file(s) that were
        selected using the file chooser dialog
    """
    if title is None:
        title = 'Select data file(s)...'
    if directory is None:
        directory = './'
    if QtGui.QApplication.instance() is None:
        app = QtGui.QApplication(sys.argv)
    fname = QtGui.QFileDialog.getOpenFileNames(None,
                                               title,
                                               directory,
                                               filter='All files (*)')
    # app.exit()

    return [str(name) for name in fname]


def shift_line_scan(s, shift, **kwargs):
    """
    shift a line scan spatially by a certain number of units

    Parameters
    ----------
    s: Signal
        s is a line scan Hyperspy signal
    shift: float
        amount by which to shift the line scan (in units of line scan,
         so usually nm)
        shift > 0 shifts a spectrum down (or right, if looking at area scan),
         while shift < 0 shifts it up
    **kwargs
        other arguments are passed to shift1D function

    Returns
    -------
    s_shifted: Signal
        shifted line scan
    """
    spectral_size = s.axes_manager[-1].size

    # create array for shift1D function
    shift_array = np.array(shift).repeat(spectral_size)

    # transpose line scan so navigation axis is signal axis
    s_transposed = s.as_spectrum(0)
    # shift by the amount in the shift_array
    s_transposed.shift1D(shift_array, **kwargs)
    # transpose back to original axes orientation
    s_shifted = s_transposed.as_spectrum(0)
    # change the title to indicate shift
    s_shifted.metadata.General.title += ' - shifted'

    return s_shifted


def shift_lines(lines,
                shifts,
                show_progressbar=True,
                progress_label='Shifting line scans:'):
    """
    Shift a list of line scans by the amount specified in shifts

    Parameters
    ----------
    lines: list of Signals
        list of Hyperspy line scans to shift
    shifts: np.array
        array (of same length as lines) that specifies shift of each line (in
        the units of the line scan)
    show_progressbar: boolean
        switch to determine whether a progressbar should be shown.
        Requires python-progressbar fork from
        https://github.com/fnoble/python-progressbar
    progress_label: str
        label to use for progressbar

    Returns
    -------
    shifted_lines: list of Signals
        duplicated list of lines, but shifted as desired
    """
    # if the length of lines and shifts do not match, raise an error
    if len(lines) is not len(shifts):
        raise ValueError("Length of lines and shifts did not match.")

    if show_progressbar:
        widgets = [progress_label, Percentage(), ' ', Bar(), ' ',
                   ETA()]
        progress = ProgressBar(widgets=widgets)

    # list to hold the shifted line scans
    shifted_lines = []

    i = 0
    try:
        prg_lines = progress(lines)
    except:
        prg_lines = lines

    # run through the lines
    for l in prg_lines:
        # shift each line by shifts[i] and add it to the shifted_lines list
        shifted_lines.append(shift_line_scan(l,
                                             shifts[i],
                                             show_progressbar=False,
                                             crop=False))
        i += 1

    return shifted_lines


def shift_area_stem(stem_s,
                    stem_list=None,
                    shifts=None,
                    crop_scan=True,
                    do_smoothing=True,
                    reset_origin=False,
                    smoothing_parameter=0.05,
                    ):
    """
    Shift a HAADF STEM signal to have a straight interface, returning a new
    Hyperspy signal with the same calibration as the original.
    Only shifts scans in the x-direction (vertical interface).
    If needed, get the shifts required from the area scan.

    Parameters
    ----------
    stem_s: HyperSpy signal with signal dimension > 1
        Should be a HAADF STEM image loaded into hyperspy that will be shifted
    stem_list: list
        list of stem signals (as output by `get_shifts_from_area_stem`.
        Providing this will save a little time extracting the scans again,
        but is honestly not very necessary
    shifts: list or None
        list of shifts to use. If None, they will be determined automatically
    crop_scan: bool
        Whether or not the resulting image should be cropped to lose the
        blank pixels resulting from shifting
    reset_origin: bool
        Whether or not to reset the origin of the x-axis to zero after cropping
    do_smoothing: bool
        Whether to smooth the signal before finding shifts
    smoothing_parameter: float
        parameter supplied to `determine_shifts` to define how much to
        smooth the data while finding the shift values
    """
    if shifts is None:
        stem_list, shifts = \
            get_shifts_from_area_stem(stem_s,
                                      do_smoothing=do_smoothing,
                                      smoothing_parameter=smoothing_parameter)

    if stem_list is None:
        stem_list = [stem_s.isig[:, i] for i in
                     range(stem_s.data.shape[0])]

    shifted_stem_list = shift_lines(stem_list, shifts)

    # Copy signal so calibration is right
    shifted_s = stem_s.deepcopy()
    shifted_s.data = hs.stack(shifted_stem_list).data
    shifted_s.metadata.General.title += ' shifted'

    if crop_scan:
        shifted_s = crop_area_scan(shifted_s, shifts)
        if reset_origin:
            shifted_s.axes_manager[0].offset = 0

    return shifted_s


def shift_area_eels(eels_s,
                    shifts,
                    crop_scan=True,
                    reset_origin=False):
    lines = [eels_s.inav[:, i] for i in range(eels_s.data.shape[0])]
    shifted_lines = shift_lines(lines,
                                shifts,
                                progress_label='Shifting EELS line scans:')

    # Copy signal so calibration is right
    shifted_s = eels_s.deepcopy()
    shifted_s.data = hs.stack(shifted_lines).data
    shifted_s.metadata.General.title += ' shifted'

    if crop_scan:
        shifted_s = crop_area_scan(shifted_s, shifts)
        if reset_origin:
            shifted_s.axes_manager[0].offset = 0

    return shifted_s



def get_shifts_from_area_stem(stem_s, **kwargs):
    """
    Get the shifts required from an area scan. Scan will be split into
    'line-scans' in the horizontal dimension

    Parameters
    ----------
    stem_s: HyperSpy signal with signal dimension > 1
        Should be a HAADF STEM image loaded into hyperspy that will be used
        to determine the needed shifts
    **kwargs:
        additional keyword arguments will be passed to `determine_shifts`

    Returns
    -------
    stem_list: list
        List of stem line scans
    shifts: np.array
        NumPy array of shift values in the x-direction (length of array is
        same as size of original signal in the y-direction)
    """
    stem_list = [stem_s.isig[:, i] for i in
                 range(stem_s.data.shape[0])]

    shifts = determine_shifts(stem_list, **kwargs)

    return stem_list, shifts


def smooth_scan(scan,
                smoothing_parm=0.05):
    """
    Smooth a signal (usually HAADF scans) and return the smoothed copy

    Parameters
    ----------
    scans: Hyperspy Signal
        scan is a Hyperspy signal that represent a STEM signal scan
    smoothing_parm: float
        amount to smooth (with Lowess filter) the scan before taking the
        derivative. The default of 0.05 seems to work well for typical STEM
        scans.

    Returns
    -------
    smoothed_scan: Hyperspy Signal
        Hyperspy signal that represent smoothed STEM signal scan
    """
    smoothed_scan = scan.deepcopy()

    smoothed_scan.smooth_lowess(smoothing_parameter=smoothing_parm,
                                number_of_iterations=1,
                                show_progressbar=False)
    return smoothed_scan


def smooth_scans(scans,
                 smoothing_parm=0.05,
                 show_progressbar=True,
                 progress_label='Smoothing scans:'):
    """
    Smooth a list of signals (usually HAADF scans) and return the smoothed copy

    Parameters
    ----------
    scans: lsit of Hyperspy Signal
        scan is list of  Hyperspy signals that represent STEM signal scans
    smoothing_parm: float or str
        amount to smooth (with Lowess filter) the scan before taking the
        derivative. The default of 0.05 seems to work well for typical STEM
        scans.
        If 'ask', a dialog will be provide to the user to get the value
    show_progressbar: boolean
        switch to determine whether or not a progress bar should be shown
        Requires python-progressbar fork from
        https://github.com/fnoble/python-progressbar
    progress_label: str
        label to use for progressbar

    Returns
    -------
    smoothed_scans: Hyperspy Signal
        Hyperspy signal that represent smoothed STEM signal scan
    """
    smoothed_scans = [0] * len(scans)

    if smoothing_parm is 'ask':
        if QtGui.QApplication.instance() is None:
            app = QtGui.QApplication(sys.argv)

        input_d = QtGui.QInputDialog(None,
                                     QtCore.Qt.WindowStaysOnTopHint)
        input_d.activateWindow()
        input_d.raise_()
        # app.exec_()

        smoothing_parm, ok = input_d.getDouble(None,
                                               'Input Dialog',
                                               'Enter desired smoothing '
                                               'parameter\n(usually '
                                               '0.03 - 0.1 works well):',
                                               decimals=2,
                                               value=0.05)
        if ok:
            # app.exit()
            pass
        else:
            print 'User cancelled input. Terminating.'
            sys.exit(1)
            # input_d.exit()

    if show_progressbar:
        widgets = [progress_label, Percentage(), ' ', Bar(), ' ', ETA()]
        progress = ProgressBar(widgets=widgets)

    i = 0
    for s in progress(scans):
        smoothed_scans[i] = smooth_scan(s, smoothing_parm)
        i += 1

    return smoothed_scans


def _interp_spectrum(stem_s, step_size, kind='cubic'):
    """
    Interpolate a signal onto a finer mesh (for subpixel accuracy)

    Parameters
    ----------
    stem_s: HyperSpy spectrum
        signal to interpolate
    step_size: float
        step size to use for interpolation
    kind: str
        kind of interpolation function to use (as defined in
        `scipy.interpolate.interp1d`)
    """
    # Create interpolation function
    f = interp1d(stem_s.axes_manager[0].axis, stem_s.data, kind=kind)

    # Define interpolation mesh
    ax_m = stem_s.axes_manager[0]
    int_range = np.arange(ax_m.axis[0], ax_m.axis[-1], step=step_size)

    # Create and fix calibration of interpolated signal
    int_s = hs.signals.Spectrum(f(int_range))
    int_s.axes_manager[0].offset = stem_s.axes_manager[0].offset
    int_s.axes_manager[0].scale = step_size
    int_s.axes_manager[0].units = stem_s.axes_manager[0].units
    int_s.axes_manager[0].name = stem_s.axes_manager[0].name
    int_s.metadata.General.title = \
        stem_s.metadata.General.title + ' - interpolated'

    return int_s


def determine_shifts(scans,
                     interp_step=0.001,
                     do_smoothing=False,
                     smoothing_parameter=0.05,
                     debug=False):
    """
    determine the units needed to shift scans to the center of spectrum image

    Parameters
    ----------
    scans: list of Signals
        scans is a list of Hyperspy signals that represent STEM signal scans
         (that match up with a series of line scans)
    interp_step: float
        step size to use for interpolation in `_interp_spectrum`
    do_smoothing: boolean
        ** Deprecated in favor of the interpolation scheme used instead **
        whether or not smoothing will be done. Disabling is useful in case
        you have already smoothed the scans ahead of time.
    smoothing_parameter: float or str
        ** Deprecated in favor of the interpolation scheme used instead **
    amount to smooth (with Lowess filter) the scan before taking the
        derivative. The default of 0.05 seems to work well for typical STEM
        scans. If 'ask', a dialog box will be risen asking for the factor
    debug: boolean
        switch whether debugging information is printed out to see the shift
        values and everything

    Returns
    -------
    shifts: np.array
        array of length len(scans) with the amount of shift (units of x-axis)
        needed to bring the transition point to the mean middle point
    """
    mids = np.zeros(len(scans))

    # smooth the scans:
    if do_smoothing:
        smoothed_scans = smooth_scans(scans,
                                      smoothing_parm=smoothing_parameter,
                                      progress_label="Smoothing STEM scans:")
        # hs.plot.plot_spectra(smoothed_scans, style='overlap')
    else:
        smoothed_scans = scans

    # Rather than derivative, let's try using another heuristic
    for i, stem in enumerate(smoothed_scans):
        interp = _interp_spectrum(stem, interp_step)

        ########################################################
        # ### New midpoint heuristic using interpolated data ###
        ########################################################
        n = 3
        first_data_i = interp.data[:n].mean()   # y-average of first N x-pixels
        last_data_i = interp.data[-n:].mean()   # y-average of last N x-pixels
        # y-average of the high and low marks
        mid_y_data_i = np.array([first_data_i, last_data_i]).mean()

        # subtract the mid y-value from the data
        subtracted_sig = interp - mid_y_data_i
        # take absolute value of this, so we can find the index of the
        # smallest difference (index of point on signal closest to the
        # midpoint)
        subtracted_sig.data = np.absolute(subtracted_sig.data)
        mid_idx = (subtracted_sig * -1).indexmax(0).data[0]
        midpoint = interp.axes_manager[0].index2value(mid_idx)

        # ### Old derivative heuristic ###
        # deriv = stem.diff(0)
        # midpoint = smoothed_scans[i].axes_manager[0].index2value(deriv
        #                                                          .indexmax(0)
        #                                                          .data)

        mids[i] = midpoint
        y_mp = interp.data[interp.axes_manager[0].value2index(
            midpoint)]

        if debug:
            # Print midpoint data and plot the data
            print("({}, {})".format(midpoint, y_mp))
            interp.plot()
            plt.scatter(midpoint, y_mp)

    shifts = mids - np.mean(mids)

    if debug:
        # print "Mids are:"
        # pprint(list(mids))
        print "Mean mids is %g" % np.mean(mids)

    # we need the negative of shifts for actual shifting
    return shifts * -1


def normalize_lines(lines,
                    show_progressbar=True,
                    progress_label='Normalizing line scans:'):
    """
    Normalize line scans by their maximum in each spectra.
    Operates in-place on the lines.

    Parameters
    ----------
    lines: list of Signals
        list of Hyperspy line scans to normalize
    show_progressbar: boolean
        switch to determine whether a progressbar should be shown.
        Requires python-progressbar fork from
        https://github.com/fnoble/python-progressbar
    progress_label: str
        label to use for progressbar
    """
    if show_progressbar:
        widgets = [progress_label, Percentage(), ' ', Bar(), ' ', ETA()]
        progress = ProgressBar(widgets=widgets)

    try:
        prg_lines = progress(lines)
    except:
        prg_lines = lines

    maxes = lines[0].max(-1).data
    try:
        # if this command is valid, we had line scans with a nav dim.
        shape = maxes.shape[0]
        line_scans = True
    except IndexError:
        # if we caught an IndexError, it means we had STEM scans that
        # did not have a navigation dimension
        line_scans = False

    for l in prg_lines:
        maxes = l.max(-1).data
        if line_scans:
            for i, s in enumerate(l):
                s.data /= maxes[i]
        else:
            for s in l:
                s.data /= maxes


def crop_area_scan(area_scan, shifts):
    """
    crop an area scan to remove areas of nan that are left from shifting lines

    Parameters:
    -----------
    area_scan: hyperspy Signal
        area spectrum image that has been built by stacking lines
    shifts: np.array
        array of shift values that were found with determine_shifts()

    Returns:
    --------
    cropped_scan: hyperspy Signal
        spatially cropped version of area_scan
    """
    # if there are more than two dimensions, we are looking at line scans
    # rather than just HAADF signals, so define the axis of interest
    # appropriately:
    if len(area_scan.data.shape) > 2:
        # here we want the x-axis, which is the first navigation axis
        axis = area_scan.axes_manager.navigation_axes[0]
    else:
        # in this case, the x-axis is actually the signal axis
        axis = area_scan.axes_manager.signal_axes[0]

    # min_shift is amount that we need to crop on the right
    # we subtract one pixel just to make sure we get rid of all
    try:
        min_shift = shifts[shifts < 0].min() - axis.scale
    except ValueError:
        min_shift = 0
    # max_shift is amount we need to crop on the left
    # again, we add one pixel to make sure we get rid of all
    max_shift = shifts[shifts > 0].max() + axis.scale

    left = max_shift
    right = axis.high_value + min_shift

    # return area_scan cropped at left and right on x-axis, and all of
    # y-axis if this is a true area spectrum image
    if len(area_scan.data.shape) > 2:
        cropped_scan = area_scan.inav[left:right, :]
    else:
        cropped_scan = area_scan.isig[left:right]

    return cropped_scan


def get_scans_and_eels_fnames():
    """
    Get the list of file names for all the line scans selected. This
    function will pop up four dialogs asking for STEM and EELS files for SiC
    to SiO2 direction and SiO2 to SiC direction. The method also sorts the
    lists in natural order so they are in an order that makes sense.

    Returns:
    --------
    c_to_o_stem, c_to_o_eels, o_to_c_stem, o_to_c_eels: list
        lists of file names for each selected files
    """
    c_to_o_stem = gui_fnames(title="Select SiC to SiO2 STEM scans...")
    c_to_o_eels = gui_fnames(title="Select SiC to SiO2 EELS lines...")
    o_to_c_stem = gui_fnames(title="Select SiO2 to SiC STEM scans...")
    o_to_c_eels = gui_fnames(title="Select SiO2 to SiC EELS lines...")

    # Sort lists to make sure they are in the right order:
    c_to_o_stem = _natural_sort(c_to_o_stem)
    c_to_o_eels = _natural_sort(c_to_o_eels)
    o_to_c_stem = _natural_sort(o_to_c_stem)
    o_to_c_eels = _natural_sort(o_to_c_eels)

    return c_to_o_stem, c_to_o_eels, o_to_c_stem, o_to_c_eels


def load_shift_and_build_area(c_to_o_stem=None,
                              c_to_o_eels=None,
                              o_to_c_stem=None,
                              o_to_c_eels=None,
                              shifts=None,
                              smoothing_parm=0.05,
                              return_unshifted=False,
                              return_uncropped=False,
                              debug=False):
    """
    Load a number of STEM signals and EELS line scans in order to
    build useful area scans out of them for decomposition and other analysis

    Usage:
    ------
    If no filenames are supplied, four file chooser dialogs will be opened.
    The files should be chosen in the order of SiC to SiO2 STEM, SiC to SiO2
    EELS, SiO2 to SiC STEM, and then SiO2 to SiC EELS.
    If there are not reversed line scans to analyze (i.e. the scans were
    acquired just in one direction), then select them in the appropriate
    window, and press 'Cancel' on the file selection for the ones that are
    not relevant.

    Note: all line scans must be same dimensions, or there will be an error.

    Parameters:
    -----------
    c_to_o_stem, c_to_o_eels, o_to_c_stem, o_to_c_eels: lists of str
        If supplied as keyword arguments, this method will not bring up a
        dialog in order to get the file names, and just use those that are
        in the lists instead. This can be useful when combined with
        ``get_scans_and_eels_fnames()`` so the function can be run multiple
        times without having to click through all the dialogs.
    shifts: list of floats
        list of shift amounts to use. Allows one to supply custom shifts for 
        each line, which will be applied to both the EELS and STEM scans
        If None, the method will try to figure it out itself
    smoothing_parameter: float or 'ask'
        This is the parameter passed to ``determine_shifts()`` in order to
        figure out how much to smooth the STEM signals before doing all the
        derivative work. Lower values are less smoothing, which will be
        more accurate, but be more susceptible to noise. Typical values are
        on the order of [0.03, 0.1], depending on the signal.
    return_unshifted: boolean
        switch whether or not to return the unshifted data (good for
        comparison)
    return_uncropped: boolean
        switch whether or not to return the uncropped data (good for
        comparison)
    debug: boolean
        switch whether debugging information is printed out to see the shift
        values and everything

    Returns:
    --------
    res: tuple
        the results tuple will have the following signals, in the following
        order:
            area_stem: hs.signal
                Hyperspy signal containing shifted and cropped STEM signals
                as an image, rather than a list of profiles
            area_eels: hs.signal
                Hyperspy signal containing the shifted and cropped EELS
                line scans as an area scan, rather than a list of single
                line scans
            file_list: list
                List of the files that were processed
            area_stem_nocrop: hs.signal
                (Optional)
                Hyperspy signal containing shifted but not cropped STEM
                signals as an image, rather than a list of profiles
            area_eels_nocrop: hs.signal
                (Optional)
                Hyperspy signal containing the shifted but not cropped EELS
                line scans as an area scan, rather than a list of single
                line scans
            area_stem_unshifted: hs.signal
                (Optional)
                Hyperspy signal containing the unshifted STEM signals as an
                image, rather than a list of profiles
            area_eels_unshifted: hs.signal
                (Optional)
                Hyperspy signal containing the unshifted EELS line scans
                as an area scan, rather than a list of single line scans
    """

    def _check_list_equal(iterator):
        # will return whether all items in list are the same or not
        return len(set(iterator)) <= 1

    # if no EELS scans are provided, get the information from dialog:
    if c_to_o_eels is None and o_to_c_eels is None:
        # get files from dialog if not supplied:
        (c_to_o_stem,
         c_to_o_eels,
         o_to_c_stem,
         o_to_c_eels) = get_scans_and_eels_fnames()

    # Save filenames in a list for reporting
    file_list = c_to_o_stem + c_to_o_eels + o_to_c_stem + o_to_c_eels

    # load in the files from the list of files:
    c_to_o_scans = [hs.load(x) for x in c_to_o_stem]
    c_to_o_lines = [hs.load(x) for x in c_to_o_eels]
    o_to_c_scans = [hs.load(x) for x in o_to_c_stem]
    o_to_c_lines = [hs.load(x) for x in o_to_c_eels]

    # flip the data in the OtoC scans and lines:
    for i in o_to_c_scans:
        i.data = i.data[::-1]
    for i in o_to_c_lines:
        i.data = i.data[::-1]

    # combine lists to make bigger lists:
    scans = c_to_o_scans + o_to_c_scans
    lines = c_to_o_lines + o_to_c_lines

    scan_sizes = [i.axes_manager.shape for i in scans]
    scan_scales = [i.axes_manager[0].scale for i in scans]
    line_sizes = [i.axes_manager.shape for i in lines]
    line_scales = [i.axes_manager[0].scale for i in lines]

    # Handle some errors related to scan sizes and magnifications
    if not _check_list_equal(scan_sizes):
        print "STEM scans were not all same size."
        print ""
        print "SiC to SiO2 files were:"
        for i in c_to_o_stem:
            print i
        print ""
        print "SiO2 to SiC files were:"
        for i in o_to_c_stem:
            print i

        print ""
        print "Sizes were:"
        pprint(scan_sizes)
        raise ValueError("All line scans must be same size for stacking.")

    if not _check_list_equal(scan_scales):
        print "STEM scans were not all same scale (different mag?)."
        print ""
        print "SiC to SiO2 files were:"
        for i in c_to_o_stem:
            print i
        print ""
        print "SiO2 to SiC files were:"
        for i in o_to_c_stem:
            print i

        print ""
        print "Scales were:"
        pprint(scan_scales)
        raise ValueError("All line scans must be same scale for stacking.")

    if not _check_list_equal(line_sizes):
        print "EELS line scans were not all same size."
        print ""
        print "SiC to SiO2 files were:"
        for i in c_to_o_eels:
            print i
        print ""
        print "SiO2 to SiC files were:"
        for i in o_to_c_eels:
            print i

        print ""
        print "Sizes were:"
        pprint(line_sizes)
        raise ValueError("All line scans must be same size for stacking.")

    if not _check_list_equal(line_scales):
        print "EELS line scans were not all same scale (different mag?)."
        print ""
        print "SiC to SiO2 files were:"
        for i in c_to_o_stem:
            print i
        print ""
        print "SiO2 to SiC files were:"
        for i in o_to_c_stem:
            print i

        print ""
        print "Scales were:"
        pprint(line_scales)
        raise ValueError("All line scans must be same scale for stacking.")

    # smooth scans:
    if shifts is None:
        smoothed_scans = smooth_scans(scans,
                                      progress_label="Smoothing STEM signals:",
                                      smoothing_parm=smoothing_parm)

    # do actual shifting and cropping:
    if shifts is None:
        shifts = determine_shifts(smoothed_scans,
                                  do_smoothing=False,
                                  debug=debug)

    if debug:
        print "Shifts are:"
        pprint(list(shifts))

    # normalize the intensity of the line scans:
    normalize_lines(lines, progress_label='Normalizing EELS line scans:')

    # normalize the intensity of the STEM profiles:
    normalize_lines(scans, progress_label='Normalizing STEM signals:')

    # shift EELS line scans
    shifted_lines = shift_lines(lines,
                                shifts,
                                progress_label='Shifting EELS line scans:')

    # shift HAADF STEM signals
    shifted_scans = shift_lines(scans,
                                shifts,
                                progress_label='Shifting STEM signals:')

    # create area spectrum images from the lines
    area_eels_nocrop = hs.utils.stack(shifted_lines)
    area_eels_nocrop.axes_manager[1].name = 'line scan'
    area_eels_nocrop.axes_manager[1].units = '#'
    area_stem_nocrop = hs.utils.stack(shifted_scans)
    area_stem_nocrop.axes_manager[0].name = 'STEM profile'
    area_stem_nocrop.axes_manager[0].units = '#'

    # Set appropriate titles for the signals
    area_eels_nocrop.metadata.General.title = 'Stacked EELS line scans - ' \
                                              'shifted'
    area_stem_nocrop.metadata.General.title = 'Stacked STEM signals - shifted'

    # crop the area spectrum images so there is no blank data
    area_eels = crop_area_scan(area_eels_nocrop, shifts)
    area_eels.axes_manager[1].name = 'line scan'
    area_eels.axes_manager[1].units = '#'
    area_stem = crop_area_scan(area_stem_nocrop, shifts)
    area_stem.axes_manager[0].name = 'STEM profile'
    area_stem.axes_manager[0].units = '#'

    # Set appropriate titles for the signals
    area_eels.metadata.General.title = 'Stacked EELS line scans - shifted ' \
                                       'and cropped'
    area_stem.metadata.General.title = 'Stacked STEM signals - shifted and ' \
                                       'cropped'

    # initialize the results list with the cropped and shifted data and the
    # list of file names that were analyzed
    res = [area_stem, area_eels, file_list]

    # if we want to return the uncropped data, add it to the list
    if return_uncropped:
        res.append(area_stem_nocrop)
        res.append(area_eels_nocrop)

    # if we want to return the unshifted data, add it to the list
    if return_unshifted:
        area_stem_unshifted = hs.utils.stack(scans)
        area_eels_unshifted = hs.utils.stack(lines)

        # Set appropriate titles for the signals
        area_eels_unshifted.metadata.General.title = 'Stacked EELS line scans'
        area_eels_unshifted.axes_manager[1].name = 'line scan'
        area_eels_unshifted.axes_manager[1].units = '#'
        area_stem_unshifted.metadata.General.title = 'Stacked STEM signals'
        area_stem_unshifted.axes_manager[0].name = 'STEM profile'
        area_stem_unshifted.axes_manager[0].units = '#'

        res.append(area_stem_unshifted)
        res.append(area_eels_unshifted)

    return res
