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


# Hyperelliptic tests over QQ
F = QQ
R.<x> = PolynomialRing(F)
Xs = [ ]

# Genus 2 factor
f = x^8 - 12*x^7 + 50*x^6 - 108*x^5 + 131*x^4 - 76*x^3 - 10*x^2 + 44*x - 19
h = R(0)
Xs.append(mHyperellipticCurve(f, h))

# Run the main functionality
for X in Xs:
    print X
    # The main functionality
    Endo = EndomorphismData(X, prec = 300, have_oldenburg = True)

    #print "Period matrix:"
    #print Endo._P_

    print "Field of definition:"
    print Endo.endomorphism_field()

    #print "Testing Rosati and degree bound:"
    #A = Endo._geo_rep_list_[1][1]
    #print A
    #print Endo.rosati_involution(A)
    #print Endo.degree_estimate(A)

    print "Geometric representation:"
    print Endo.geometric().representation()

    print "Over several fields:"
    #print Endo.geometric().representation()
    #print Endo.over_base().representation()
    R.<t> = PolynomialRing(F)
    K.<s> = NumberField(t^2 - 2)
    overK = Endo.over_field(K)
    print K
    #print overK.representation()
    #print overK.algebra()
    #print overK.description()
    print overK.pretty_print()

    print "Lattice:"
    #print Endo.lattice()
    #print Endo.lattice().representations()
    #print Endo.lattice().algebras()
    #print Endo.lattice().descriptions()
    print Endo.lattice().pretty_print()

    #print "Verification:"
    #print Endo.dimension_algebra()
    #print Endo.base_point()
    #A = Endo._geo_rep_dict_[2]['tangent']
    #print A
    #print Endo.correspondence(A)
    #print Endo.verify()

    #print "Testing same functionality over a field:"
    #A = overK._list_[1][1][1]
    #print A
    #print overK.rosati_involution(A)
    #print overK.degree_estimate(A)
    #print overK.dimension_algebra()
    #print overK.verify_algebra()
    #print overK.verify_saturated()
    #print overK.base_point()
    #print overK.correspondence(A)
    #print overK.verify()
    #print overK.full()

    #print "Decomposition:"
    Dec = Endo.decomposition()
    print Dec
    print Dec.field
    print Dec.idempotents()
    #print Dec.projections()
    print Dec.factors()
    #print Dec._factors_desc_()
    #print Dec.verify()
