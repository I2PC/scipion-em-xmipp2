# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Carlos Oscar Sanchez Sorzano
# *              Estrella Fernandez Gimenez
# *
# *  BCU, Centro Nacional de Biotecnologia, CSIC
# *
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

import os
from pyworkflow.em.convert import ImageHandler
import pyworkflow.utils as pwutils
from pyworkflow.utils.path import makePath

from tomo.protocols.protocol_base import ProtTomoSubtomogramAveraging
from pyworkflow.protocol.params import PointerParam
from xmipp2.convert import writeSetOfVolumes


class Xmipp2ProtMLTomo(ProtTomoSubtomogramAveraging):
    """ Protocol to align subtomograms using MLTomo """

    _label = 'mltomo'
        
    def __init__(self, **args):
        ProtTomoSubtomogramAveraging.__init__(self, **args)


    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input subtomograms')
        form.addParam('inputVolumes', PointerParam, pointerClass="SetOfSubTomograms", label='Set of volumes',
                      help="Set of subtomograms to align with MLTomo")
        form.addParallelSection(threads=0, mpi=4)

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        #self._insertFunctionStep('runMLTomo')
        #self._insertFunctionStep('createOutput')

    #--------------------------- STEPS functions -------------------------------
    def convertInputStep(self):
        fnDir=self._getExtraPath("inputVolumes")
        makePath(fnDir)
        fnRoot=os.path.join(fnDir,"subtomo")
        writeSetOfVolumes(self.inputVolumes.get(),fnRoot)
        self.runJob("xmipp_selfile_create",'"%s*.vol">%s'%(fnRoot,self._getExtraPath("subtomograms.sel")))

    def runMLTomo(self):
        outputFile = self._getExtraPath(pwutils.replaceBaseExt(inputMic, "txt"))  # Cambiar????
        args = " --image=%s --outfile=%s" % (inputMic, outputFile)
        self.runJob("xmipp_ml_tomo", args, numberOfMPIs=self.numberOfMpi.get())
    
    def createOutput(self):
        pass

    #--------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        return summary


    def _methods(self):
        methods = []
        methods.append(" ")
        return methods
