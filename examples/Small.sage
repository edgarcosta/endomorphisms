"""
 *  Some examples of endomorphism rings
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

load("../Initialize.sage")

F = QQ
R.<x> = PolynomialRing(F)
f, h = [x^6+3*x^5+2*x^4+7*x^3+11*x^2+14, x^2+x]
f, h = [x^7+x^6-4*x^5+4*x^3-5*x^2+2*x-1, x^4+x^3+x+1]
f, h = [-3*x^7-3*x^6-3*x^5+x^2+x, x^3+1]
#f, h = [-x^8+3*x^7-3*x^6-2*x^5-2*x^4+x^2, x^4+x^3+x+1]
#f, h = [x^5 -3*x^4 -2*x - 1, x^3 + x^2 + x + 1]
#f, h = [x^4 + 7, x^3 + x]
#f, h = [x^4 - 7, x^3 + x]
#f, h = [x^4 + x^2, x^3 + 1]
#f, h = [ -38*x^6 + 102*x^5 + 77*x^4 - 276*x^3 - 118*x^2 + 232*x + 132, 0 ]

X = mHyperellipticCurve(f, h)

print X
Endo = EndomorphismData(X, prec = 300, have_oldenburg = True)

print ""
print "Field of definition:"
print Endo.endomorphism_field()

print ""
print "Lattice:"
print Endo.lattice().pretty_print()

Dec = Endo.decomposition()
print ""
print "Decomposition:"
print Dec.field
print Dec.factors()
print Dec.verify()
