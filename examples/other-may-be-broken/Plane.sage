"""
 *  Plane examples (need extra code to be run)
"""

from endomorphisms import EndomorphismData, mHyperellipticCurve, mPlaneCurve

F = QQ
R.<x,y,z> = PolynomialRing(F)
f = y^3*z - x^4 - z^4
f = x^4 + 2*x^2*y^2 + 2*x^2*y*z - 2*x^2*z^2 + 2*y^4 + 4*y^3*z - y^2*z^2 - 3*y*z^3 + z^4
X = mPlaneCurve(f, 200)

Endo = EndomorphismData(X)
lat = Endo.lattice()
print lat.pretty_print()

#print ""
#print lat._desc_
#print lat._sthash_

#Dec = Endo.decomposition()
#print ""
#print Dec.factors()
