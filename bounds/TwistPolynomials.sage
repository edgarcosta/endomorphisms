"""
 *  Bound functionality
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

def twistPolynomial(poly, k) :
    R.<x,v> = PolynomialRing(ZZ, 2)
    f = R(poly)
    g = v-x^k
    T.<v> = PolynomialRing(QQ)
    return T(f.resultant(g))
