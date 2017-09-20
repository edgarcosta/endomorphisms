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
f, h = [x^5-x^4-2*x^3+x^2+x,0]
f, h = [-6*x^6-4*x^5-4*x^4-5*x^3-2*x^2-x-1,1]
f, h = [x^5-3*x^4+3*x^2+x,0]
f, h = [7*x^4+72*x^2+219,x^3+x]
#f, h = [10*x^6-24*x^5+34*x^4-29*x^3+17*x^2-6*x+1,1]
X = mHyperellipticCurve(f, h)

#R.<x,y,z> = PolynomialRing(F)
#f = x^3*z + x*y^3 + x*z^3 + y^4 + y^3*z
#f = x^3*z + x^2*z^2 + x*y^2*z + x*y*z^2 + x*z^3 - y^4 - 2*y^3*z - 2*y^2*z^2 - y*z^3
#f = x^3*z + x^2*y*z + x^2*z^2 + x*y^3 + x*y^2*z + x*y*z^2 + x*z^3 + y^4 + y^3*z
#X = mPlaneCurve(f)

print ""
print "Curve:"
print X
Endo = EndomorphismData(X, prec = 300, have_oldenburg = True)

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

print ""
print "Endomorphisms:"
print Endo.lattice().pretty_print()
#print "Verify:"
#test = Endo.verify()
#print test
