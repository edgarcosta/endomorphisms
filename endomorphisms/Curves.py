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

from sage.all import magma

def mHyperellipticCurve(f, h):
    return magma.HyperellipticCurve(magma(f), magma(h))

def mPlaneCurve(f):
    return magma.PlaneCurve(magma(f))
