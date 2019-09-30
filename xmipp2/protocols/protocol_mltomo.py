# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez
# *              Carlos Oscar Sanchez Sorzano
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
from os.path import exists

from pyworkflow.em.data import Transform
from pyworkflow.utils.path import makePath
from pyworkflow.protocol.params import PointerParam, BooleanParam, IntParam, StringParam, LEVEL_ADVANCED

from tomo.objects import AverageSubTomogram
from tomo.protocols.protocol_base import ProtTomoSubtomogramAveraging
from xmipp2.convert import writeVolume, writeSetOfVolumes, eulerAngles2matrix


class Xmipp2ProtMLTomo(ProtTomoSubtomogramAveraging):
    """ Protocol to align subtomograms using MLTomo. It only supports alignment on the axis Y.
        MLTomo aligns and classifies 3D images with missing data regions in Fourier
        space, e.g. subtomograms or RCT reconstructions, by a 3D
        multi-reference refinement based on a maximum-likelihood (ML) target
        function."""

    _label = 'mltomo'

    def __init__(self, **args):
        ProtTomoSubtomogramAveraging.__init__(self, **args)


    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input subtomograms')
        form.addParam('inputVolumes', PointerParam, pointerClass="SetOfSubTomograms",
                      label='Set of volumes', help="Set of subtomograms to align with MLTomo")
        form.addParam('randomInitialization', BooleanParam, default=True,
                      label='Random initialization of classes:', help="Initialize randomly the first classes. If you "
                           "don't initialize randomly, you must supply a set of initial classes")
        form.addParam('initialClasses', PointerParam, label="Initial classes",
                      condition="not randomInitialization", pointerClass='SetOfClassesSubTomograms',
                      help='Set of initial classes to start the classification')
        form.addParam('numberOfReferences', IntParam, label='Number of references', default=10,
                      condition="randomInitialization", help="Number of references to generate automatically")
        form.addParam('numberOfIters', IntParam, label='Number of iterations', default=15,
                      help="Number of iterations to perform")
        form.addParam('downscDim', IntParam, label='Downscaled dimension', default=30,
                      help="Use downscaled (in fourier space) images of this size")
        form.addParam('angularSampling', IntParam, label='Angular sampling rate', default=15,
                      help="Angular sampling rate (in degrees)")
        form.addParam('inputMask', PointerParam, label="Mask",
                      allowsNull=True, pointerClass='VolumeMask',
                      help='Optionally, select a mask. If a mask is used, the program will keep '
                      'rotations and translations from docfile fixed, only classify')
        form.addParam('extraParams', StringParam, label='Extra parameters', default='-perturb -dont_impute',
                      expertLevel=LEVEL_ADVANCED,
                      help="Any of the parameters of the MLTomo (https://github.com/I2PC/xmipp-portal/wiki/Ml_tomo )")
        form.addParallelSection(threads=0, mpi=8)

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runMLTomo')
        self._insertFunctionStep('createOutput')

    #--------------------------- STEPS functions -------------------------------
    def convertInputStep(self):
        fnDir=self._getExtraPath("inputVolumes")
        makePath(fnDir)
        fnRoot=os.path.join(fnDir,"subtomo")
        writeSetOfVolumes(self.inputVolumes.get(),fnRoot)
        self.fnSel= self._getExtraPath("subtomograms.sel")
        self.runJob("xmipp_selfile_create",'"%s*.vol">%s'%(fnRoot,self.fnSel),numberOfMpi=1)
        if self.initialClasses.get() is not None:
            fnRootRef=os.path.join(fnDir,"reference")
            writeSetOfVolumes(self.initialClasses.get().iterRepresentatives(),fnRootRef)
            self.fnSelRef= self._getExtraPath("references.sel")
            self.runJob("xmipp_selfile_create",'"%s*.vol">%s'%(fnRootRef,self.fnSelRef),numberOfMpi=1)
        if self.inputMask.get() is not None:
            self.fnMask = os.path.join(fnDir, "mask.vol")
            writeVolume(self.inputMask.get(), self.fnMask)

    def runMLTomo(self):
        self._createFilesForMLTomo()
        args = ' -i ' + self._getExtraPath("subtomograms.sel") + \
               ' -o ' + self._getExtraPath("mltomo") + \
               ' -doc ' + self._getExtraPath("subtomograms.doc") + \
               ' -iter ' + str(self.numberOfIters.get()) + \
               ' -ang ' + str(self.angularSampling.get()) + \
               ' -dim ' + str(self.downscDim.get()) + \
               ' ' + self.extraParams.get()
        if self.initialClasses.get() is not None:
            args = args + ' -ref ' + self.fnSelRef
        else:
            args = args + ' -nref ' + str(self.numberOfReferences.get())
        if self.inputMask.get() is not None:
            args = args + ' -mask ' + self.fnMask + ' -dont_align'
        fhWedge = self._getExtraPath("wedge.doc")
        if exists(fhWedge):
            args = args + ' -missing ' + fhWedge
        self.runJob("xmipp_ml_tomo", args, numberOfMpi=self.numberOfMpi.get())
    
    def createOutput(self):
        self.subtomoSet = self._createSetOfSubTomograms()
        inputSet = self.inputVolumes.get()
        self.subtomoSet.copyInfo(inputSet)
        if self.numberOfIters < 10:
            self.fnDoc = '%s/mltomo_it00000%d.doc' % (self._getExtraPath(),self.numberOfIters)
        else:
            self.fnDoc = '%s/mltomo_it0000%d.doc' % (self._getExtraPath(),self.numberOfIters)
        self.docFile = open(self.fnDoc)
        self.subtomoSet.copyItems(inputSet, updateItemCallback=self._updateItem)
        self.docFile.close()
        classesSubtomoSet = self._createSetOfClassesSubTomograms(self.subtomoSet)
        classesSubtomoSet.classifyItems(updateClassCallback=self._updateClass)
        self._defineOutputs(outputSubtomograms=self.subtomoSet)
        self._defineSourceRelation(self.inputVolumes, self.subtomoSet)
        self._defineOutputs(outputClassesSubtomo=classesSubtomoSet)
        self._defineSourceRelation(self.inputVolumes, classesSubtomoSet)
        # self._cleanFiles()

    #--------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputClassesSubtomo'):
            summary.append("Input subtomograms: *%d* \nRequested classes: *%d*\nGenerated classes: *%d* in *%d* iterations\n"
                       % (self.inputVolumes.get().getSize(),self.numberOfReferences,
                          self.outputClassesSubtomo.getSize(),self.numberOfIters))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputClassesSubtomo'):
            methods.append('We classified %d subtomograms from %s into %d classes %s using *MLTomo*.'
                       % (self.inputVolumes.get().getSize(),self.getObjectTag('inputVolumes'),
                          self.outputClassesSubtomo.getSize(),self.getObjectTag('outputClassesSubtomo')))
        else:
            methods.append("Output classes not ready yet.")
        return methods

    def _citations(self):
        return ['Scheres2009c']

    #--------------------------- UTILS functions ----------------------------------
    def _createFilesForMLTomo(self):
        mw = 0
        if (self.inputVolumes.get().getFirstItem().getAcquisition().getAngleMin()):
            mw = 1
            wedgeDict={}
            for subtomogram in self.inputVolumes.get():
                key=(subtomogram.getAcquisition().getAngleMin(),subtomogram.getAcquisition().getAngleMax())
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

    def _updateItem(self, item, row):
        nline = self.docFile.next()
        if nline.startswith(' ;'):
            nline = self.docFile.next()
        if nline.startswith(' ;'):
            nline = self.docFile.next()
        nline = nline.rstrip()
        id = int(nline.split()[0])
        if (item.getObjId() == id):
            rot = nline.split()[2]
            tilt = nline.split()[3]
            psi = nline.split()[4]
            shiftx = nline.split()[5]
            if shiftx.startswith('-'):
                shiftx = shiftx[1:]
            else:
                shiftx = '-' + shiftx
            shifty = nline.split()[6]
            if shifty.startswith('-'):
                shifty = shifty[1:]
            else:
                shifty = '-' + shifty
            shiftz = nline.split()[7]
            if shiftz.startswith('-'):
                shiftz = shiftz[1:]
            else:
                shiftz = '-' + shiftz
            refId = float(nline.split()[8])
            A = eulerAngles2matrix(rot, tilt, psi, shiftx, shifty, shiftz)
            transform = Transform()
            transform.setMatrix(A)
            item.setTransform(transform)
            item.setClassId(refId)

    def _updateClass(self, item):
        classId = item.getObjId()
        item.setAlignment3D()
        directory = self._getExtraPath()
        fnRep = ('%s/mltomo_ref00000%d.vol' % (directory, classId))
        representative = AverageSubTomogram()
        representative.setLocation(1, fnRep)
        representative.copyInfo(self.subtomoSet)
        representative.setClassId(classId)
        item.setRepresentative(representative)

    def _cleanFiles(self):
        for iter in range(0, int(self.numberOfIters)+1):
            os.remove(self._getExtraPath('mltomo_it00000%d.sel' % iter))
            for ref in range(1, int(self.numberOfReferences)+1):
                os.remove(self._getExtraPath('mltomo_it00000%d_ref00000%d.vol' % (iter, ref)))
                os.remove(self._getExtraPath('mltomo_it00000%d_wedge00000%d.vol' % (iter, ref)))
        for iter in range(1, int(self.numberOfIters)+1):
            os.remove(self._getExtraPath('mltomo_it00000%d.fsc' % iter))

        for ref in range(1, int(self.numberOfReferences) + 1):
            if os.stat(self._getExtraPath('mltomo_it00000%d_ref00000%d.sel' % (self.numberOfIters, ref))).st_size == 0:
                os.remove(self._getExtraPath('mltomo_it00000%d_ref00000%d.sel' % (self.numberOfIters, ref)))
                os.remove(self._getExtraPath('mltomo_ref00000%d.vol' % (ref)))
                continue
            else:
                continue
        for iter in range(1, int(self.numberOfIters)):
            os.remove(self._getExtraPath('mltomo_it00000%d.doc' % iter))
            for ref in range(1, int(self.numberOfReferences) + 1):
                os.remove(self._getExtraPath('mltomo_it00000%d_ref00000%d.sel' % (iter,ref)))
