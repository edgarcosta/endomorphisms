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

def CertifyNonIsogenous(LPolys1, LPolys2, geometric = true) :
    extDegree = 12
    if not geometric :
        extDegree = 1
    for p in range (2,maxP):
      if is_prime(p) and LPolys1[p] <> 0 and LPolys2[p] <> 0 :
        q1 = LPolys1[p]
        q2 = LPolys2[p]

        # print "Twisting the L-polys"
        q1 = twistPolynomial(q1, extDegree)
        q2 = twistPolynomial(q2, extDegree)

        if q1 != q2 :
            return true

    return false
