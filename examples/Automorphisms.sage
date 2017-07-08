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

R.<x> = PolynomialRing(QQ)
f = x^3 - 1
h = 0
X = mHyperellipticCurve(f, h)

prec = 300
CC = magma.ComplexField(prec)
s = magma.Sqrt(CC(-21))

# Examples by Christophe Ritzenthaler
Ps = [ ]

P = magma.Transpose(magma.Matrix([
[1, 0,   5,  2*s],
[0, 1, -2*s,   17]
]))
Ps.append(P)

P = magma.Transpose(magma.Matrix([
[1, 0,        25,  12*s - 5],
[0, 1, -12*s - 5,       122]
]))
Ps.append(P)

P = magma.Transpose(magma.Matrix([
[1, 0,         25,  12*s - 10],
[0, 1, -12*s - 10,        125]
]))
Ps.append(P)

P = magma.Transpose(magma.Matrix([
[1, 0,        125,  62*s + 20],
[0, 1, -62*s + 20,        649]
]))
Ps.append(P)


for P in Ps:
    Endo = EndomorphismData(X, prec = prec, have_oldenburg = False, periods = P)

    print "Field of definition:"
    print Endo.endomorphism_field()

    print "Geometric representation:"
    Geo = Endo.geometric()
    rep = Geo.representation()
    print Geo.pretty_print()
    Rsfixed = Geo.rosati_fixed_module()
    print magma.MinimalPolynomial(Rsfixed[2])
