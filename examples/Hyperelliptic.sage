"""
 *  Hyperelliptic examples
"""

from endomorphisms import EndomorphismData, mHyperellipticCurve, mPlaneCurve

F = QQ
R.<x> = PolynomialRing(F)
f, h = [21*x^7 + 37506*x^5 + 933261*x^3 + 5841759*x, R(0)]
f, h = [x^7 + 6*x^5 + 9*x^3 + x, R(0)]
f, h = [16*x^7 + 357*x^5 - 819*x^3 + 448*x, R(0)]
f, h = [-4*x^8 + 105*x^6 - 945*x^4 + 2100*x^2 - 5895*x + 420, x^4]
f, h = [x^7 - 14*x^6 + 210*x^5 - 658*x^4 + 245*x^3 + 588*x^2 + 637*x - 686, R(0)]
# Conductor 2^4 13^4, all twists have rational point, CM curve
f, h = [x^6 - 8*x^4 - 8*x^3 + 8*x^2 + 12*x - 8, R(0)]
X = mHyperellipticCurve(f, h, 300)

Endo = EndomorphismData(X)
lat = Endo.lattice()
print lat.pretty_print()

# TODO: This needs an update
#Dec = Endo.decomposition()
#print ""
#print Dec.factors()
