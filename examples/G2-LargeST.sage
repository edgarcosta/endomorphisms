"""
 *  An example with the largest possible Sato-Tate group in genus 2
"""

from endomorphisms import EndomorphismData, mHyperellipticCurve

F = QQ
R.<x> = PolynomialRing(F)
Xs = [ ]

# Largest genus 2 case (takes a LOT of time)
f = x^6 - 5*x^4 + 10*x^3 - 5*x^2 + 2*x - 1
h = R(0)
Xs.append(mHyperellipticCurve(f, h))

for X in Xs:
    print X
    Endo = EndomorphismData(X, prec = 300, molin_neurohr = False)

    print "Field of definition:"
    print Endo.endomorphism_field()

    print "Geometric representation:"
    print Endo.geometric().representation()

    print "Over several fields:"
    R.<t> = PolynomialRing(F)
    K.<s> = NumberField(t^2 - 2)
    overK = Endo.over_field(K)
    print K
    print overK.pretty_print()

    print "Lattice:"
    print Endo.lattice().pretty_print()

    Dec = Endo.decomposition()
    print Dec
    print Dec.field
    print Dec.idempotents()
    print Dec.factors()
