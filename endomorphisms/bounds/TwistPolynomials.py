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

from sage.all import ZZ, QQ, PolynomialRing

def twistPolynomial(poly, k) :
    R = PolynomialRing(ZZ, 'x,v')
    x, v =  R.gens()
    f = R(poly)
    g = v-x^k
    T = PolynomialRing(QQ, 'v')
    return T(f.resultant(g))
