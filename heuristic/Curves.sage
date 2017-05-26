"""
 *  Creation of curves
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

def mHyperellipticCurve(f, h):
    return magma.HyperellipticCurve(magma(f), magma(h))

def mPlaneCurve(f):
    mf = magma(f)
    mP2 = magma.ProjectiveSpace(magma.Parent(mf))
    return magma.Curve(magma.Scheme(mP2, mf))
