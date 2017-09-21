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
f, h = [x^4 + x^2, x^3 + 1]
X = mHyperellipticCurve(f, h)

R.<x,y,z> = PolynomialRing(F)
f = x^4 + 2*x^3*z - x^2*y^2 + 2*x^2*y*z - x^2*z^2 - x*y^3 - 2*x*y^2*z + x*y*z^2 - 2*x*z^3 + y^4 - y^3*z + 2*y^2*z^2 - y*z^3 + z^4
X = mPlaneCurve(f)

print ""
print "Curve:"
print X
Endo = EndomorphismData(X, prec = 100, have_oldenburg = True)

Dec = Endo.decomposition()
print ""
print "Decomposition:"
Dec = Endo.decomposition()
print "Factors:"
facs = Dec.factors()
print facs
print "Verify:"
test = Dec.verify()
print test
print Dec.correspondences()

print ""
print "Endomorphisms:"
print Endo.lattice().pretty_print()
#print "Verify:"
#test = Endo.verify()
#print test
