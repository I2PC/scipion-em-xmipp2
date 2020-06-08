# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez
# *              Carlos Oscar Sanchez Sorzano
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
This module implements visualization program
for ml_tomo.
"""
import os
from xmipp2.protocols import Xmipp2ProtMLTomo
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from pyworkflow.protocol.params import LabelParam
from pwem.viewers.plotter import EmPlotter
from pwem.viewers.views import DataView


class Xmipp2ProtMlTomoViewer(ProtocolViewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj. """

    _label = 'viewer Ml_Tomo'
    _targets = [Xmipp2ProtMLTomo]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def _defineParams(self, form):
        form.addSection(label='Resolution')
        form.addParam('doShowFsc', LabelParam,
                      label="Display Fourier Shell Correlation")

    def _getVisualizeDict(self):
        return {'doShowFsc': self._viewFsc}

    def _viewFsc(self, e=None):
        return self._loadPlots("Fourier Shell Correlation", 'FSC', color='r')

    def _loadPlots(self, title, plotLabel, **kwargs):
        """ Check if the FSC metadata is generated and if so,
        read the plots and the metadata.
        *args and **kwargs will be passed to self._createPlot function.
        """
        fnFsc = self.protocol._getExtraPath("mltomo.fsc")

        if not os.path.exists(fnFsc):
            return [self.errorMessage('The FSC file was not produced\n',
                                      title='Missing result file')]
        # leer archivo FSC
        fnFsc = open(self.protocol._getExtraPath("mltomo.fsc"),'r')
        lines = fnFsc.readlines()
        freq = []
        # fnClasses = self.protocol._getExtraPath("Classes") # generated classes
        classes = [] # j,i (j num classes, i num lines)
        for i in lines:
            freq.append(i.split()) # column frecs
            # classes.append(i.split()[j]) # j num class

        fnFsc.close()

        plotter = EmPlotter(self, freq, classes[1], title, **kwargs) # change class index
        return [plotter, DataView(plotLabel)]