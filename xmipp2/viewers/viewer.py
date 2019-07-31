# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from pyworkflow.em.viewers import showj
from pyworkflow.viewer import Viewer, DESKTOP_TKINTER, WEB_DJANGO, CommandView
from pyworkflow.em.protocol import *

from pyworkflow.em.viewers.viewers_data import DataViewer
from pyworkflow.em.viewers.plotter import EmPlotter
from pyworkflow.em.viewers.views import CtfView, ObjectView
from pyworkflow.em.viewers.showj import *
from pyworkflow.em.viewers.viewer_monitors import MovieGainMonitorPlotter

import xmippLib
from xmipp2.convert import *
from xmipp2.protocols.protocol_mltomo import Xmipp2ProtMLTomo
import numpy as np



class XmippViewer(DataViewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [
        Xmipp2ProtMLTomo
    ]

    def __createTemporaryCtfs(self, obj, setOfMics):
        pwutils.cleanPath(obj._getPath("ctfs_temporary.sqlite"))
        self.protocol._createFilenameTemplates()
        ctfSet = self.protocol._createSetOfCTF("_temporary")

        for mic in setOfMics:
            micDir = obj._getExtraPath(removeBaseExt(mic.getFileName()))
            ctfparam = self.protocol._getFileName('ctfparam', micDir=micDir)

            if exists(ctfparam) or exists('xmipp_default_ctf.ctfparam'):
                if not os.path.exists(ctfparam):
                    ctfparam = 'xmipp_default_ctf.ctfparam'
                # ctfModel = readCTFModel(ctfparam, mic)
                # self.protocol._setPsdFiles(ctfModel, micDir)
                # ctfSet.append(ctfModel)

        if not ctfSet.isEmpty():
            ctfSet.write()
            ctfSet.close()

        return ctfSet

    def _visualize(self, obj, **kwargs):
        cls = type(obj)
        if issubclass(cls, SetOfNormalModes):
            fn = obj.getFileName()
            # from .viewer_nma import OBJCMD_NMA_PLOTDIST, OBJCMD_NMA_VMD
            # objCommands = "'%s' '%s'" % (OBJCMD_NMA_PLOTDIST, OBJCMD_NMA_VMD)
            # self._views.append(ObjectView(self._project, self.protocol.strId(),
            #                               fn, obj.strId(),
            #                               viewParams={OBJCMDS: objCommands},
            #                               **kwargs))


        else:
            # Use default visualization defined in base class
            DataViewer._visualize(self, obj, **kwargs)

        return self._views

    def getCTFViews(self, ctfSet):
        # This could be used by any CTF viewer to show CTF plus, phaseShift plot
        # if applies.
        # Return phaseShift plot if apply
        firstCtf = ctfSet.getFirstItem()

        if firstCtf.hasPhaseShift():
            phase_shift = []

            for ctf in ctfSet.iterItems():
                phShift = ctf.getPhaseShift()
                phase_shift.append(phShift)

            plotter = EmPlotter()
            plotter.createSubPlot("Phase Shift estimation",
                                  "Number of CTFs", "Phase Shift")
            plotter.plotData(np.arange(0, len(phase_shift)), phase_shift)
            self._views.append(plotter)

        # Return Standard CTF view (showJ)
        self._views.append(CtfView(self._project, ctfSet))
