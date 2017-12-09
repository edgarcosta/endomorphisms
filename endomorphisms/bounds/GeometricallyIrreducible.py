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

def IsGeometricallyIrreducible(LPolys):
    for p in range (2,maxP):                                # loop over primes
      if LPolys[p] <> 0 :               # ensure p is of good reduction
       # print "Testing p =",p
       q = LPolys[p]
       g = q.degree() / 2
       q = twistPolynomial(q, extensionBounds[g])
       if(q.coefficients(sparse=False)[g] % p<>0) :         # check for ordinarity
            if(q.is_irreducible()) :
                return True;
    return False;
