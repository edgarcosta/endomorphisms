"""
 *  An example from the paper
"""

from endomorphisms import EndomorphismData, mHyperellipticCurve

# Ambient ring:
F = QQ
CCSmall = magma.ComplexField(5)
R.<x> = PolynomialRing(F)

# Curve:
f = x^5 - x^4 + 4*x^3 - 8*x^2 + 5*x - 1
h = R(0)
X = mHyperellipticCurve(f, h)

print X
Endo = EndomorphismData(X, prec = 100, have_oldenburg = False)

print "Period matrix:"
P = magma.Transpose(Endo._P_)
print magma.ChangeRing(P, CCSmall)

print "Endomorphism over QQ (sqrt (2)):"
R.<t> = PolynomialRing(F)
K.<s> = NumberField(t^2 - 2)
overK = Endo.over_field(K)
endodict = overK.full()
M = endodict['representation'][1]['tangent']
MCC = endodict['representation'][1]['approx']
R = endodict['representation'][1]['homology']
print M
print R

#CC = magma.BaseRing(P)
#MCC = magma.ChangeRing(MCC, CC)
#RCC = magma.ChangeRing(R, CC)
#print magma.ChangeRing(MCC * P - P * RCC, CCSmall)

print "Verifying endomorphisms:"
test = Endo.verify()
print test
