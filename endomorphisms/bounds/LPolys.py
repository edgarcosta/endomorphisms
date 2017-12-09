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

'''
Compute the characteristic polynomial of Frobenius for the curve y^2=f(x) at the prime p.
It assumes that p is odd and the curve has good reduction at p
'''

from sage.all import PolynomialRing, ZZ

def Frobenius_Polynomial(p,C) :
    R = PolynomialRing(GF(p), 't')
    pCount=C.count_points_exhaustive(n=2)
    a=pCount[0]
    b=pCount[1]
    s1=p+1-a;
    s2=1/2*(a^2+2*p-2*a-2*a*p+b)
    poly=ZZ['x']([p^2,-p*s1,s2,-s1,1])
    return poly
