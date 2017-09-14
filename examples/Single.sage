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

R.<x,y,z> = PolynomialRing(F)
f = x^3*z + x^2*y*z + x^2*z^2 + x*y^3 - 3*x*y^2*z - 4*x*z^3 - y^4 + 2*y^3*z + 2*z^4
X = mPlaneCurve(f)

print X
Endo = EndomorphismData(X, prec = 100, have_oldenburg = True)
print Endo.verify()

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
