"""
 *  Some examples of genus 3 curves with RM
"""

from endomorphisms import EndomorphismData, mHyperellipticCurve

F = QQ
R.<x> = PolynomialRing(F)

f = x^7 - 14*x^6 + 210*x^5 - 658*x^4 + 245*x^3 + 588*x^2 + 637*x - 686
h = 0
X = mHyperellipticCurve(f, h)

print X
Endo = EndomorphismData(X, prec = 300, molin_neurohr = False)

print "Field of definition:"
print Endo.endomorphism_field()

print "Geometric representation:"
overK = Endo.geometric()
rep = overK.representation()
print overK.has_generator(B = 10)
print overK.few_generators()
