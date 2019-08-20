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
import fnmatch
from os.path import exists

from pyworkflow.em.data import Transform
from pyworkflow.utils.path import makePath
from pyworkflow.protocol.params import PointerParam, IntParam, StringParam, LEVEL_ADVANCED

from tomo.objects import SubTomogram, AverageSubTomogram
from tomo.protocols.protocol_base import ProtTomoSubtomogramAveraging
from xmipp2.convert import writeSetOfVolumes, eulerAngles2matrix


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
        self._insertFunctionStep('createOutput')

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
        self._cleanFiles()
        directory = self._getExtraPath()
        self._classesInfo = []
        subtomoSet = self._createSetOfSubTomograms()
        subtomoSet.setSamplingRate(self.inputVolumes.get().getSamplingRate())
        self.referencesSet = self._createSetOfAverages()
        self.referencesSet.setSamplingRate(self.inputVolumes.get().getSamplingRate())
        fnDoc = '%s/mltomo_it00000%d.doc' % (directory,self.numberOfIters)
        for refId in range(0,self.numberOfReferences):
            for filename in os.listdir(directory):
                if fnmatch.fnmatch(filename, 'mltomo_it00000%d_ref00000%d.sel' % (self.numberOfIters,refId)):
                    name = '%s/mltomo_it00000%d_ref00000%d.sel' % (directory,self.numberOfIters,refId)
                    if not os.stat(name).st_size == 0:
                        f = open(name)
                        reference = AverageSubTomogram()
                        reference.setFileName('%s/mltomo_ref00000%d.vol' % (directory, refId))
                        reference.setClassId(refId)
                        reference.setSamplingRate(self.inputVolumes.get().getSamplingRate())
                        self.referencesSet.append(reference)
                        if not refId in self._classesInfo:
                            self._classesInfo.append(refId)
                        for line in f:
                            line = line.rstrip()
                            fnSubtomo = line.split(None, 1)[0]
                            subtomo = SubTomogram()
                            subtomo.setFileName(fnSubtomo)
                            subtomo.setSamplingRate(self.inputVolumes.get().getSamplingRate())
                            subtomo.setClassId(refId)
                            transform = Transform()
                            d = open(fnDoc)
                            for dline in d:
                                if fnSubtomo in dline:
                                    nline = d.next()
                                    nline = nline.rstrip()
                                    rot = nline.split()[2]
                                    tilt = nline.split()[3]
                                    psi = nline.split()[4]
                                    shiftx = nline.split()[5]
                                    if shiftx.startswith('-'): shiftx = shiftx[1:]
                                    else: shiftx = '-' + shiftx
                                    shifty = nline.split()[6]
                                    if shifty.startswith('-'): shifty = shifty[1:]
                                    else: shifty = '-' + shifty
                                    shiftz = nline.split()[7]
                                    if shiftz.startswith('-'): shiftz = shiftz[1:]
                                    else: shiftz = '-' + shiftz
                            d.close()
                            A = eulerAngles2matrix(rot, tilt, psi, shiftx, shifty, shiftz)
                            transform.setMatrix(A)
                            subtomo.setTransform(transform)
                            subtomoSet.append(subtomo)
                            if not line:
                                continue
                        f.close()
                    continue
                else:
                    continue

        classesSubtomoSet = self._createSetOfClassesSubTomograms(subtomoSet)
        classesSubtomoSet.classifyItems(updateClassCallback=self._updateClass)
        self._defineOutputs(outputSubtomograms=subtomoSet)
        self._defineSourceRelation(self.inputVolumes, subtomoSet)
        self._defineOutputs(outputClassesSubtomo=classesSubtomoSet)
        self._defineSourceRelation(self.inputVolumes, classesSubtomoSet)

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

    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            item.setAlignment3D()
            directory = self._getExtraPath()
            fn = ('%s/mltomo_ref00000%d.vol' % (directory, classId))
            item.getRepresentative().setLocation(classId,fn)

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
