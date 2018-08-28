"""
 *  Plane examples (need extra code to be run)
"""

from endomorphisms import EndomorphismData, mHyperellipticCurve, mPlaneCurve

F = QQ
R.<x> = PolynomialRing(F)
# Conductor 2^4 13^4, all twists have rational point, CM curve
f, h = [x^6 - 8*x^4 - 8*x^3 + 8*x^2 + 12*x - 8, 0];
X = mHyperellipticCurve(f, h)

Endo = EndomorphismData(X, prec = 200, molin_neurohr = True)
lat = Endo.lattice()
stdesc = lat._stdesc_
print lat.pretty_print()

Dec = Endo.decomposition()
print Dec.factors()
