# **************************************************************************
# *
# * Authors:    Estrella Fernandez Gimenez [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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

from pyworkflow.tests import BaseTest, setupTestProject
from tomo.protocols import ProtImportSubTomograms
from tomo.tests import DataSet
from xmipp2.protocols import Xmipp2ProtMLTomo, Xmipp2ProtKerdensom


class TestXmipp2Mltomo(BaseTest):
    """This class check if the protocol Mltomo in xmipp2 works properly"""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.setOfSubtomograms = cls.dataset.getFile('basename.hdf')

    def _runMltomo(self, randomInitialization=True, numberOfReferences=2,
                   numberOfIters=3, angularSampling=15):
        protImport = self.newProtocol(ProtImportSubTomograms,
                                      filesPath=self.setOfSubtomograms,
                                      samplingRate=5)
        self.launchProtocol(protImport)

        protMltomo = self.newProtocol(Xmipp2ProtMLTomo,
                                      inputVolumes=protImport.outputSubTomograms,
                                      randomInitialization=randomInitialization,
                                      numberOfReferences=numberOfReferences,
                                      numberOfIters=numberOfIters,
                                      angularSampling=angularSampling)
        self.launchProtocol(protMltomo)
        self.assertIsNotNone(protMltomo.outputSubtomograms,
                             "There was a problem with SetOfSubtomograms output")
        self.assertIsNotNone(protMltomo.outputClassesSubtomo,
                             "There was a problem with outputClassesSubtomo output")
        return protMltomo

    def test_outputMltomo(self):
        protMltomo = self._runMltomo()
        outputSubtomos = getattr(protMltomo, 'outputSubtomograms')
        outputClasses = getattr(protMltomo, 'outputClassesSubtomo')
        self.assertTrue(outputSubtomos)
        self.assertTrue(outputSubtomos.getFirstItem().hasTransform())
        self.assertTrue(outputClasses)
        self.assertTrue(outputClasses.hasRepresentatives())
        return protMltomo

class TestXmippKerdensom(BaseTest):
    """This class check if the protocol to calculate the kerdensom from particles in Xmipp works properly."""
    @classmethod
    def setData(cls, dataProject='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.particlesFn = cls.dataset.getFile('particles')
        cls.particlesDir = cls.dataset.getFile('particlesDir')
        cls.volumesFn = cls.dataset.getFile('volumes')
        cls.volumesDir = cls.dataset.getFile('volumesDir')
        cls.averagesFn = cls.dataset.getFile('averages')
        cls.averagesDir = cls.dataset.getFile('averagesDir')

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        BaseTest.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.align2D = cls.runCL2DAlign(cls.protImport.outputParticles)

    def test_kerdensom(self):
        print("Run Kerdensom")
        xmippProtKerdensom = self.newProtocol(Xmipp2ProtKerdensom, SomXdim=2, SomYdim=2)
        xmippProtKerdensom.inputParticles.set(self.align2D.outputParticles)
        self.launchProtocol(xmippProtKerdensom)
        self.assertIsNotNone(xmippProtKerdensom.outputClasses, "There was a problem with Kerdensom")

    def test_kerdensomMask(self):
        print("Run Kerdensom with a mask")
        protMask = self.runCreateMask(3.5, 100)
        xmippProtKerdensom = self.newProtocol(Xmipp2ProtKerdensom, SomXdim=2, SomYdim=2)
        xmippProtKerdensom.inputParticles.set(self.align2D.outputParticles)
        xmippProtKerdensom.useMask.set(True)
        xmippProtKerdensom.Mask.set(protMask.outputMask)
        self.launchProtocol(xmippProtKerdensom)
        self.assertIsNotNone(xmippProtKerdensom.outputClasses, "There was a problem with Kerdensom")