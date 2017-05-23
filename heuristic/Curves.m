/***
 *  Properties of curves
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic CurveType(X::Crv) -> MonStgElt
{String that gives type.}

if Type(X) eq CrvHyp then
    return "hyperelliptic";
elif Type(X) eq CrvPln then
    return "plane";
end if;

end intrinsic;


intrinsic EmbedCurveEquations(X::Crv, prec::RngIntElt) -> MonStgElt
{Gives equations of a curve as elements of CC.}

if CurveType(X) eq "hyperelliptic" then
    f, h := HyperellipticPolynomials(X);
    return EmbedAsComplexPolynomials([ f, h ], prec);
elif CurveType(X) eq "plane" then
    F := DefiningPolynomial(X);
    return EmbedAsComplexPolynomials([ F ], prec);
end if;

end intrinsic;
