# **************************************************************************
# *
# * Authors:  Laura del Cano (ldelcano@cnb.csic.es)
# *           Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es) [1]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
This package contains the protocols and data for APPION
"""
import os
import pyworkflow.em

from pyworkflow.utils import Environ
from .constants import DOGPICKER_HOME, V0_2_1_1


_references = ['Voss2009']
_logo = 'appion_logo.png'


class Plugin(pyworkflow.em.Plugin):
    _homeVar = DOGPICKER_HOME
    _pathVars = [DOGPICKER_HOME]
    _supportedVersions = V0_2_1_1

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(DOGPICKER_HOME, 'dogpicker-0.2.1.1')

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch Appion. """
        environ = Environ(os.environ)

        environ.update({
            'PATH': Plugin.getHome(),
            'LD_LIBRARY_PATH': str.join(cls.getHome(), 'appionlib')
                               + ":" + cls.getHome(),
        }, position=Environ.BEGIN)

        return environ

    @classmethod
    def isVersionActive(cls):
        return cls.getActiveVersion().startswith(V0_2_1_1)

    @classmethod
    def defineBinaries(cls, env):

        # Add dogpicker
        env.addPackage('dogpicker', version='0.2.1.1',
                       tar='dogpicker-0.2.1.1.tgz',
                       default=True)


pyworkflow.em.Domain.registerPlugin(__name__)





