"""
 *  Plane examples (need extra code to be run)
"""

from endomorphisms import EndomorphismData, mHyperellipticCurve, mPlaneCurve

R.<x,y,z> = PolynomialRing(F)
f = y^3*z - x^4 - z^4
X = mPlaneCurve(f)

Endo = EndomorphismData(X, prec = 200, molin_neurohr = True)
lat = Endo.lattice()
stdesc = lat._stdesc_
print lat.pretty_print()

Dec = Endo.decomposition()
print Dec.factors()
