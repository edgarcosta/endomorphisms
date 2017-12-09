"""
 *  An example from the paper
"""

from endomorphisms import EndomorphismData, mHyperellipticCurve

# Ambient ring:
F = QQ
R.<x> = PolynomialRing(F)

# Curve:
f = -3*x^6 + 8*x^5 - 30*x^4 + 50*x^3 - 71*x^2 + 50*x - 27
h = R(0)
X = mHyperellipticCurve(f, h)

print X
Endo = EndomorphismData(X, prec = 300, have_oldenburg = False)

print "Endomorphism over QQ (sqrt (5)):"
R.<t> = PolynomialRing(F)
K.<s> = NumberField(t^2 - 5)
overK = Endo.over_field(K)
endodict = overK.full()
print endodict['representation'][1]['tangent']
print endodict['representation'][1]['homology']

print "Test for saturation:"
print Endo.verify_saturated()

print "Verifying endomorphisms:"
test = Endo.verify()
print test
