"""
 *  Relative splitting fields with GP optimizations
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

from sage.all import magma, PolynomialRing, QQ

def Relative_Splitting_Field_Extra(fs, bound = 0):
    bound_set = (bound != 0)
    F = magma.BaseRing(fs[1])
    overQQ = (magma.Degree(F) == 1)
    if overQQ:
        R = PolynomialRing(QQ, 'x')
    fs = sorted(fs, key = lambda f : -magma.Degree(f))
    K = F
    for f in fs:
        if not magma.HasRoot(f, K):
            for tup in magma.Factorization(f, K):
                K = magma.ExtendRelativeSplittingField(K, F, tup[1], Optimize = False)
                if overQQ:
                    g = magma.DefiningPolynomial(K)
                    g = R(str(gp.polredabs(R(magma.Eltseq(g)))))
                    K = magma.NumberField(g)
                else:
                    K = magma.ClearFieldDenominator(K)
                if bound_set and magma.Degree(K) >= bound:
                    K = magma.MakeExtension(K, F)
                    K = magma.DefineOrExtendInfinitePlaceFunction(K);
                    return K
    K = magma.MakeExtension(K, F)
    K = magma.DefineOrExtendInfinitePlaceFunction(K);
    return K
