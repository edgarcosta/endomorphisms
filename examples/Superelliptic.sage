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


# Plane tests
F = QQ
P2.<x,y,z> = ProjectiveSpace(F, 2)

R.<t> = PolynomialRing(F)
p = t^4 + 2*t^2 + 3*t + 5

f = y^4 - (x^4 + 2*x^2*z^2 + 3*x*z^3 + 5*z^4)
f = y^3*z - (x^4 + 2*x^2*z^2 + 3*x*z^3 + 5*z^4)
f = y^3*z - (x^4 + 2*x^2*z^2 + 5*z^4)
X = mPlaneCurve(f)

prec = 300

print X
Endo = EndomorphismData(X, prec = prec, have_oldenburg = True)

#print "Period matrices:"
#print "Directly:"
#print magma.Transpose(magma.SE_BigPeriodMatrix(magma(p), 3, Prec = prec))
#print "Via amended functionality:"
#print Endo._P_

print "Field of definition:"
print Endo.endomorphism_field()

print "Geometric representation:"
print Endo.geometric().representation()
