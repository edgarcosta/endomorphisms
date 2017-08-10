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

# Plane test
F = QQ
P2.<x,y,z> = ProjectiveSpace(F, 2)
Xs = [ ]

f = x^4 - x^3*y + 2*x^3*z + 2*x^2*y*z + 2*x^2*z^2 - 2*x*y^2*z + 4*x*y*z^2 - y^3*z + 3*y^2*z^2 + 2*y*z^3 + z^4
Xs.append(mPlaneCurve(f))

for X in Xs:
    print X
    Endo = EndomorphismData(X, prec = 50, have_oldenburg = True)

    print "Field of definition:"
    print Endo.endomorphism_field()

    print "Geometric representation:"
    print Endo.geometric().representation()

    print "Lattice:"
    print Endo.lattice().pretty_print()
    Endo.lattice().optimize_representations()
    print Endo.lattice().representations()
