"""
 *  An example from the paper
"""

from endomorphisms import EndomorphismData, mHyperellipticCurve

# Ambient ring:
F = QQ
R.<x> = PolynomialRing(F)

# Curve:
f = x^8 - 12*x^7 + 50*x^6 - 108*x^5 + 131*x^4 - 76*x^3 - 10*x^2 + 44*x - 19
h = R(0)
X = mHyperellipticCurve(f, h)

print X
Endo = EndomorphismData(X, prec = 300, molin_neurohr = False)

print "Geometric representation:"
overK = Endo.over_field("geometric")
endodict = overK.full()
rep = overK.representation()
print rep

print "Verifying endomorphisms:"
print Endo.verify()

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
