# **************************************************************************
# *
# * Authors:  Estrella Fernandez Gimenez
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
1. Write from base classes to Xmipp2.4 specific files
2. Read from Xmipp2.4 files to base classes
"""
import numpy as np
from pyworkflow.em.convert import ImageHandler


def writeSetOfVolumes(setOfVolumes, outputFnRoot):
    ih = ImageHandler()
    i = 1
    for volume in setOfVolumes:
        ih.convert(volume, "%s%06d.vol"%(outputFnRoot,i))
        i+=1

def eulerAngles2matrix(alpha, beta, gamma, shiftx, shifty, shiftz):
    A = np.empty([4,4])
    A.fill(2)
    A[3,3] = 1
    A[3,0:3] = 0

    A[0,3] = float(shiftx)
    A[1,3] = float(shifty)
    A[2,3] = float(shiftz)

    alpha = float(alpha)
    beta = float(beta)
    gamma = float(gamma)

    sa = np.sin(np.deg2rad(alpha))
    ca = np.cos(np.deg2rad(alpha))
    sb = np.sin(np.deg2rad(beta))
    cb = np.cos(np.deg2rad(beta))
    sg = np.sin(np.deg2rad(gamma))
    cg = np.cos(np.deg2rad(gamma))

    cc = cb * ca
    cs = cb * sa
    sc = sb * ca
    ss = sb * sa

    A[0,0] = cg * cc - sg * sa
    A[0,1] = cg * cs + sg * ca
    A[0,2] = -cg * sb
    A[1,0] = -sg * cc - cg * sa
    A[1,1] = -sg * cs + cg * ca
    A[1,2] = sg * sb
    A[2,0] = sc
    A[2,1] = ss
    A[2,2] = cb

    return A
