"""
 *  Plane examples (need extra code to be run)
"""

from endomorphisms import EndomorphismData, mHyperellipticCurve, mPlaneCurve

F = QQ
R.<x> = PolynomialRing(F)
f, h = [x^4 + x^2, x^3 + 1]
#f, h = [x^5 - x^4 + 4*x^3 - 8*x^2 + 5*x - 1, 0]
X = mHyperellipticCurve(f, h)

R.<x,y,z> = PolynomialRing(F)
f = x^4 + 2*x^3*z - x^2*y^2 + 2*x^2*y*z - x^2*z^2 - x*y^3 - 2*x*y^2*z + x*y*z^2 - 2*x*z^3 + y^4 - y^3*z + 2*y^2*z^2 - y*z^3 + z^4
f = x^3*z + x^2*y*z + x^2*z^2 - x*y^3 + x*y^2*z + x*z^3 - y^2*z^2 + y*z^3
f = x^3*z + x^2*y^2 + x^2*y*z + x*y^3 + x*y^2*z + x*y*z^2 + x*z^3 + y^3*z + y^2*z^2
X = mPlaneCurve(f)

print ""
print "Curve:"
print X
Endo = EndomorphismData(X, prec = 300, molin_neurohr = True)

print ""
print "Endomorphisms:"
print Endo.lattice().pretty_print()
