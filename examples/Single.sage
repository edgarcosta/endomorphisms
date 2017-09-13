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
f, h = x^6 + 3*x^5 + 2*x^4 + 7*x^3 + 11*x^2 + 14, x^2 + x
X = mHyperellipticCurve(f, h)

#R.<x,y,z> = PolynomialRing(F)
#f = x^3*z + x^2*y*z + x^2*z^2 - x*y^3 - x*y*z^2 - x*z^3 + y^2*z^2
#f = x^3*z + y^4 + y^3*z + y^2*z^2 + y*z^3
#X = mPlaneCurve(f)

print X
Endo = EndomorphismData(X, prec = 100, have_oldenburg = True)

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
print Dec.idempotents()
print Dec.verify()
