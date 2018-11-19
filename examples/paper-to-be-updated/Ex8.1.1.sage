"""
 *  An example from the paper
"""

from endomorphisms import EndomorphismData, mHyperellipticCurve

# Ambient ring:
F = QQ
R.<x> = PolynomialRing(F)

# Curve:
f = 5*x^6 + 10*x^3 - 4*x + 1
h = R(0)
X = mHyperellipticCurve(f, h)

print X
Endo = EndomorphismData(X, prec = 300, molin_neurohr = False)

print "Endomorphism over QQ (sqrt (5)):"
R.<t> = PolynomialRing(F)
K.<s> = NumberField(t^2 - 5)
overK = Endo.over_field(K)
endodict = overK.full()
print endodict['representation'][0]['tangent'] + endodict['representation'][1]['tangent']
print endodict['representation'][0]['homology'] + endodict['representation'][1]['homology']

print "Verifying endomorphisms:"
test = Endo.verify()
print test
