"""
 *  An example from the paper
"""

from endomorphisms import EndomorphismData, mHyperellipticCurve

# Ambient ring:
F = QQ
R.<x> = PolynomialRing(F)

# Curve:
f = 24*x^5 + 36*x^4 - 4*x^3 - 12*x^2 + 1
h = R(0)
X = mHyperellipticCurve(f, h)

print X
Endo = EndomorphismData(X, prec = 300, have_oldenburg = False)

print "Field of definition:"
print Endo.endomorphism_field()

print "Endomorphism over QQ (sqrt (-3)):"
R.<t> = PolynomialRing(F)
K.<s> = NumberField(t^2 + 3)
overK = Endo.over_field(K)
endodict = overK.full()
print endodict['representation'][1]['tangent']
print endodict['representation'][1]['homology']
print "Degree estimate:"
print Endo.degree_estimate(endodict['representation'][1]['tangent'])

print "Verifying endomorphisms:"
test = Endo.verify()
print test
