"""
 *  An example from the paper
"""

from endomorphisms import EndomorphismData, mHyperellipticCurve

# Ambient ring:
F = QQ
R.<x> = PolynomialRing(F)

# Curve:
f = x^6 - 8*x^4 + 2*x^3 + 16*x^2 - 36*x - 55
h = R(0)
X = mHyperellipticCurve(f, h)

print X
Endo = EndomorphismData(X, prec = 300, have_oldenburg = False)

print "Decomposition:"
Dec = Endo.decomposition()
print Dec.field

print Dec.idempotents()
print Dec.projections()
print Dec.factors()
print Dec.verify()
print Dec.correspondences()

print "Verifying endomorphisms:"
test = Endo.verify()
print test
