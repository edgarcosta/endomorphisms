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

def certifySurfaceNonQM(LPolys, conductor):
    g = 2                                                   # only for surfaces
    for p in range (2,maxP):                                # loop over primes
      if is_prime(p) and conductor % p <> 0 :               # ensure p is of good reduction
        q = LPolys(p)
        q = twistPolynomial(q, extensionBounds[g])

        if(q.is_square() == false) :
            return true
    return false;
