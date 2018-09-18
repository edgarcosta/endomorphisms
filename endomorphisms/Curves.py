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

def mHyperellipticCurve(f, h, prec):
    return magma.HyperellipticCurveExtra(magma(f), magma(h), magma(prec))

def mPlaneCurve(f, prec):
    return magma.PlaneCurveExtra(magma(f), magma(prec))
