# **************************************************************************
# *
# * Authors:  Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
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
import math
import numpy as np
from pwem.objects import Transform
from pwem.emlib.image import ImageHandler


def writeVolume(volume, outputFn):
    ih = ImageHandler()
    ih.convert(volume, "%s" % outputFn)

def writeSetOfVolumes(setOfVolumes, outputFnRoot):
    ih = ImageHandler()
    for volume in setOfVolumes:
        i = volume.getObjId()
        ih.convert(volume, "%s%06d.vol"%(outputFnRoot,i))

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

def matrix2eulerAngles(A):
    abs_sb = np.sqrt(A[0, 2] * A[0, 2] + A[1, 2] * A[1, 2])
    if (abs_sb > 16*np.exp(-5)):
        gamma = math.atan2(A[1, 2], -A[0, 2])
        alpha = math.atan2(A[2, 1], A[2, 0])
        if (abs(np.sin(gamma)) < np.exp(-5)):
            sign_sb = np.sign(-A[0, 2] / np.cos(gamma))
        else:
            if np.sin(gamma) > 0:
                sign_sb = np.sign(A[1, 2])
            else:
                sign_sb = -np.sign(A[1, 2])
        beta = math.atan2(sign_sb * abs_sb, A[2, 2])
    else:
        if (np.sign(A[2, 2]) > 0):
            alpha = 0
            beta  = 0
            gamma = math.atan2(-A[1, 0], A[0, 0])
        else:
            alpha = 0
            beta  = np.pi
            gamma = math.atan2(A[1, 0], -A[0, 0])
    gamma = np.rad2deg(gamma)
    beta  = np.rad2deg(beta)
    alpha = np.rad2deg(alpha)
    return alpha, beta, gamma, A[0,3], A[1,3], A[2,3]


def readDocfile(self, item):
    nline = self.docFile.readline()
    if nline.startswith(' ;'):
        nline = self.docFile.readline()
    if nline.startswith(' ;'):
        nline = self.docFile.readline()
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

def writeDocfile(self, fhSel, fhDoc, volumes, wedge):
    fhDoc.write(" ; Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Zoff (6), Ref (7), Wedge (8), "
                "Pmax/sumP (9), LL (10)\n")
    for line in fhSel:
        imgName = line.split()[0]
        for vol in volumes.iterItems():
            if ('%06d' % vol.getObjId()) in line:
                rot, tilt, psi, xoff, yoff, zoff = matrix2eulerAngles(vol.getTransform().getMatrix())
                xoff = xoff * (-1)
                yoff = yoff * (-1)
                zoff = zoff * (-1)
                classid = vol.getClassId()
                if classid is None:
                    classid = 0
                fhDoc.write(" ; %s\n%d 10 %f %f %f %f %f %f %d %d 0 0\n" % (imgName, vol.getObjId(), rot, tilt, psi,
                                                                            xoff, yoff, zoff, classid, wedge))
