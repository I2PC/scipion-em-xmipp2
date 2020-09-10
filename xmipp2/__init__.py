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
This package contains the protocols and data for xmipp 2.4
"""
import os
import pwem
from pyworkflow.utils import Environ
from .constants import XMIPP2_HOME

__version__ = '3.0.1'

class Plugin(pwem.Plugin):
    _homeVar = XMIPP2_HOME
    _pathVars = [XMIPP2_HOME]


    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(XMIPP2_HOME, "Xmipp-2.4-src")


    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch Xmipp2. """
        environ = Environ(os.environ)

        environ.update({
            'PATH': os.path.join(cls.getVar(XMIPP2_HOME),'bin'),
            'LD_LIBRARY_PATH': os.path.join(cls.getVar(XMIPP2_HOME),'lib'),
        }, position=Environ.BEGIN)

        return environ

    @classmethod
    def isVersionActive(cls):
        return cls.getActiveVersion().startswith("")

    @classmethod
    def defineBinaries(cls, env):
        compileCmd = [("./scons.configure && ./scons.compile -j 8",["lib/libXmippRecons_Interface.so","bin/xmipp_mpi_ml_tomo"])]
        env.addPackage('xmipp2', version='2.4',
                       tar='Xmipp-2.4-src.tgz',
                       commands=compileCmd,
                       default=True)


# @classmethod
# def defineBinaries(cls, env):
#     """ Define the Xmipp binaries/source available tgz.
#     """
#     scons = tryAddPipModule(env, 'scons', '3.0.4', default=True)
#
# def tryAddPipModule(env, moduleName, *args, **kwargs):
#     """ To try to add certain pipModule.
#         If it fails due to it is already add by other plugin or Scipion,
#           just returns its name to use it as a dependency.
#         Raise the exception if unknown error is gotten.
#     """
#     try:
#         return env.addPipModule(moduleName, *args, **kwargs)._name
#     except Exception as e:
#         if "Duplicated target '%s'" % moduleName == str(e):
#             return moduleName
#         else:
#             raise Exception(e)