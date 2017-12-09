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

from sage.all import NumberField

def IsGeometricallyNotSquareOfCM(LPolys):
    g = 2						    # only for surfaces
    CMDisc = -1
    for p in range (2,maxP):                                # loop over primes
      if LPolys[p] <> 0 :
        q = LPolys[p]
        q = twistPolynomial(q, 2*extensionBounds[g])
        if(q.is_square() == False) :
            return True
	r = q.factor()
	r = r[0][0];
	if r.degree() == 2 :
		if r.coefficients(sparse = False)[1] % p <> 0 :
			K = NumberField(r, 'u');
			d = K.disc()
		        if CMDisc == -1 :
				CMDisc = d
			if CMDisc <> d :
				return True
    return False;
