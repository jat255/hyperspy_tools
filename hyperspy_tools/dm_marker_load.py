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

from hyperspy import api as hs
import hyperspy.io_plugins.digital_micrograph as dm


def plot_survey_with_markers(fname, add_text=True):
    """
    Plot a hyperspy signal with the markers from digital micrograph enabled

    Parameters
    ----------
    fname : str
        Name of .dm3 file to load
    add_text : str
        Switch to control if labels for the markers are added to the plot
    """

    def _add_beam(image, location):
        beam_m = hs.plot.markers.point(x=location[1],
                                       y=location[0],
                                       color='red')
        image.add_marker(beam_m)
        if add_text:
            beam_text_m = hs.plot.markers.text(x=location[1] - 0.5,
                                               y=location[0] - 1.5,
                                               color='red',
                                               text='Beam',
                                               size='xx-small')
            image.add_marker(beam_text_m)

    def _add_si(image, location):
        # adds a green rectangle (or line, if the coordinates are such) to
        # image
        si_m = hs.plot.markers.rectangle(x1=location[1],
                                         y1=location[0],
                                         x2=location[3],
                                         y2=location[2],
                                         color='#13FF00')
        image.add_marker(si_m)
        if add_text:
            si_text_m = hs.plot.markers.text(x=location[1],
                                             y=location[0] - 0.5,
                                             color='#13FF00',
                                             text='Spectrum Image',
                                             size='xx-small')
            image.add_marker(si_text_m)

    def _add_drift(image, location):
        drift_m = hs.plot.markers.rectangle(x1=location[1],
                                            y1=location[0],
                                            x2=location[3],
                                            y2=location[2],
                                            color='yellow')
        image.add_marker(drift_m)
        if add_text:
            drift_text_m = hs.plot.markers.text(x=location[1],
                                                y=location[0] - 0.5,
                                                color='yellow',
                                                text='Spatial Drift',
                                                size='xx-small')
            image.add_marker(drift_text_m)

    im = hs.load(fname)
    flist = dm.file_reader(fname)
    annotation_list = flist[0]['original_metadata']['DocumentObjectList'][
                               'TagGroup0']['AnnotationGroupList']

    mapping = {
        'Beam': _add_beam,
        'Spectrum Image': _add_si,
        'Spatial Drift': _add_drift
    }

    im.plot(cmap='gist_gray')

    for i in range(len(annotation_list)):
        try:
            label = annotation_list['TagGroup' + str(i)]['Label']
            loc = annotation_list['TagGroup' + str(i)]['Rectangle']
            scaled_loc = [im.axes_manager[0].scale * i for i in loc]
            mapping[label](im, scaled_loc)

        except KeyError:
            pass
