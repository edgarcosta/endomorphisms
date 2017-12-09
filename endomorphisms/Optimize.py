"""
 *  Improve representation of endomorphisms
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

from sage.all import magma

def Optimize_Representation(rep):
    K = magma.BaseRing(rep[1])
    F = magma.BaseRing(K)
    if (not magma.IsQQ(F)) or (magma.Degree(K) == 1):
        return rep
    R = magma.PolynomialRing(magma.Rationals())
    g = magma.DefiningPolynomial(K)
    g = R(str(gp.polredabs(R(magma.Eltseq(g)))))
    L = magma.NumberField(g)
    return magma.TransferMatrices(rep, K, L)
