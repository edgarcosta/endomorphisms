"""
 *  Representation functions
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

def repr_curve(X):
    curve_type = magma.CurveType(X)
    if str(curve_type) == "hyperelliptic":
        f, h = magma.HyperellipticPolynomials(X, nvals = 2)
        if magma.IsZero(h):
            return " the hyperelliptic curve y^2 = {}".format(str(f))
        else:
            return " the hyperelliptic curve y^2 + ({})*y = {}".format(str(h), str(f))
    elif str(curve_type) == "plane":
        F = magma.DefiningPolynomial(X)
        return " the plane curve {} = 0".format(str(F))

def repr_endomorphism_data(End):
    return "The endomorphism data of" + repr_curve(End.X)

def repr_lattice(Lat):
    return "The endomorphism lattice of" + repr_curve(Lat.X)

def repr_over_field(over_field):
    pre = "The endomorphism structure of" + repr_curve(over_field.X)
    if over_field.field == "geometric":
        post = " over the algebraic closure of its base field"
    elif over_field.field == "base":
        post = " over its base field"
    else:
        post = " over " + str(over_field.field)
    return pre + post

def repr_decomposition(decomp):
    return "The decomposition structure of" + repr_curve(decomp.X)
