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
from os.path import exists

from tomo.protocols.protocol_base import ProtTomoSubtomogramAveraging
from pyworkflow.protocol.params import PointerParam, IntParam, StringParam, LEVEL_ADVANCED
from xmipp2.convert import writeSetOfVolumes



class Xmipp2ProtMLTomo(ProtTomoSubtomogramAveraging):
    """ Protocol to align subtomograms using MLTomo. It only supports alignment on the axis Y """

    _label = 'mltomo'

    def __init__(self, **args):
        ProtTomoSubtomogramAveraging.__init__(self, **args)


    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input subtomograms')
        form.addParam('inputVolumes', PointerParam, pointerClass="SetOfSubTomograms",
                      label='Set of volumes', help="Set of subtomograms to align with MLTomo")
        form.addParam('numberOfReferences', IntParam, label='Number of references', default=10,
                      help="Number of references to generate automatically")
        form.addParam('numberOfIters', IntParam, label='Number of iterations', default=15,
                      help="Number of iterations to perform")
        form.addParam('angularSampling', IntParam, label='Angular sampling rate', default=15,
                      help="Angular sampling rate (in degrees)")
        form.addParam('downscDim', IntParam, label='Downscaled dimension', default=30,
                      help="Use downscaled (in fourier space) images of this size")
        form.addParam('extraParams', StringParam, label='Extra parameters', default='-perturb -dont_impute',
                      expertLevel=LEVEL_ADVANCED,
                      help="Any of the parameters of the MLTomo (https://github.com/I2PC/xmipp-portal/wiki/Ml_tomo )")
        form.addParallelSection(threads=0, mpi=8)

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runMLTomo')
        #self._insertFunctionStep('createOutput')

    #--------------------------- STEPS functions -------------------------------
    def convertInputStep(self):
        fnDir=self._getExtraPath("inputVolumes")
        makePath(fnDir)
        fnRoot=os.path.join(fnDir,"subtomo")
        writeSetOfVolumes(self.inputVolumes.get(),fnRoot)
        self.fnSel= self._getExtraPath("subtomograms.sel")
        self.runJob("xmipp_selfile_create",'"%s*.vol">%s'%(fnRoot,self.fnSel),numberOfMpi=1)

    def createFilesForMLTomo(self):
        mw = 0
        if (self.inputVolumes.get().getFirstItem().getAcquisitionAngleMin()):
            mw = 1
            wedgeDict={}
            for subtomogram in self.inputVolumes.get():
                key=(subtomogram.getAcquisitionAngleMin(),subtomogram.getAcquisitionAngleMax())
                if not key in wedgeDict:
                    wedgeDict[key]=[]
                wedgeDict[key].append(subtomogram)

            fhWedge = open(self._getExtraPath("wedge.doc"),'w')
            fhWedge.write(" ; Wedgeinfo\n ; wedge_y\n")
            i=1
            for key in wedgeDict:
                fhWedge.write("%d 2 %d %d\n" %(i,key[0],key[1]))
                i+=1
            fhWedge.close()

        fhDoc = open(self._getExtraPath("subtomograms.doc"),'w')
        j = 1
        with open(self._getExtraPath("subtomograms.sel"),'r') as fhSel:
            for line in fhSel:
                imgName = line.split()[0]
                fhDoc.write(" ; %s\n%d 10  0 0 0 0 0 0 0 %d 0 0\n" %(imgName,j,mw))
                j+=1
        fhDoc.close()
        fhSel.close()

    def runMLTomo(self):
        fnIn = self._getExtraPath("subtomograms.sel")
        fnOut = self._getExtraPath("mltomo")
        self.createFilesForMLTomo()
        fhDoc = self._getExtraPath("subtomograms.doc")
        args = ' -i ' + fnIn + \
               ' -o ' + fnOut + \
               ' -nref ' + str(self.numberOfReferences.get()) + \
               ' -doc ' + fhDoc + \
               ' -iter ' + str(self.numberOfIters.get()) + \
               ' -ang ' + str(self.angularSampling.get()) + \
               ' -dim ' + str(self.downscDim.get()) + \
               ' ' + self.extraParams.get()

        fhWedge = self._getExtraPath("wedge.doc")
        if exists(fhWedge):
            args = args + ' -missing ' + fhWedge

        self.runJob("xmipp_ml_tomo", args, numberOfMpi=self.numberOfMpi.get())
    
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
