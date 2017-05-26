"""
 *  Some examples of endomorphism rings
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

load("../Initialize.sage")

Xs = [ ]

# RM over QQ
#R.<x> = PolynomialRing(QQ)
#f,h = [2*x^5-3*x^4+x^3+x^2-x,1]
#Xs.append(mHyperellipticCurve(f, h))

# Example by Haluk
R.<t> = PolynomialRing(QQ)
F.<u> = NumberField(t^2 - t + 1)
R.<x> = PolynomialRing(F)
f = R([ -30*u + 42, -156*u + 312, -66*u + 186, -1456*u + 1040, -90*u + 126, 156*u - 312, -22*u + 62 ])
h = R(0)
Xs.append(mHyperellipticCurve(f, h))

for X in Xs:
    print X
    # The main functionality
    Endo = EndomorphismData(X, prec = 300, have_oldenburg = True)

    print "Period matrix:"
    print Endo.period_matrix()

    print "Field of definition:"
    print Endo.endomorphism_field()

    print "Geometric representation:"
    Geo = Endo.geometric()
    print Geo.representation()
    print Geo.algebra()
    print Geo.pretty_print()
