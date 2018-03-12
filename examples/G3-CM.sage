"""
 *  Some examples of genus 3 curves with CM
"""

from endomorphisms import EndomorphismData, mHyperellipticCurve

F = QQ
R.<x> = PolynomialRing(F)

f = 21*x^7 + 37506*x^5 + 933261*x^3 + 5841759*x
h = 0
f = x^7 + 6*x^5 + 9*x^3 + x
h = 0
f = 16*x^7 + 357*x^5 - 819*x^3 + 448*x
h = 0
f = -4*x^8 + 105*x^6 - 945*x^4 + 2100*x^2 - 5895*x + 420
h = x^4
X = mHyperellipticCurve(f, h)

print X
Endo = EndomorphismData(X, prec = 300, molin_neurohr = False)

print "Field of definition:"
print Endo.endomorphism_field()

print "Geometric representation:"
overK = Endo.geometric()
print overK.representation()
print overK.has_generator()
print overK.few_generators()
print overK.pretty_print()
