# **************************************************************************
# *
# * Authors:  Laura del Cano (ldelcano@cnb.csic.es)
# *           Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es)
# *
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
This module contains converter functions that will serve to:
1. Write from base classes to Appion specific files
2. Read from Appion files to base classes
"""

from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.data import Coordinate
from pyworkflow.em.metadata import MetaData, MDL_XCOOR, MDL_YCOOR

from pyworkflow.utils.path import replaceBaseExt, join, exists
import xmippLib

def readCoordinates(mic, fileName, coordsSet):
    if exists(fileName):
         md = MetaData()
         md.readPlain(fileName, 'xcoor ycoor')
         for objId in md:
            x = md.getValue(MDL_XCOOR, objId)
            y = md.getValue(MDL_YCOOR, objId)
            coord = Coordinate()
            coord.setPosition(x, y)
            coord.setMicrograph(mic)
            coordsSet.append(coord)


def writeSetOfCoordinates():
    pass


def writeSetOfVolumes(setOfVolumes, outputFnRoot):
    ih = ImageHandler()
    i = 1
    for volume in setOfVolumes:
        ih.convert(volume, "%s%06d.vol"%(outputFnRoot,i))
        i+=1
























