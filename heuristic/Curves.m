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


intrinsic PlaneCurve(f::RngMPolElt) -> Crv
{Returns the plane curve defined by f.}

R := Parent(f); gens := GeneratorsSequence(R);
if #gens eq 3 then
    return Curve(Scheme(ProjectiveSpace(R), f));
end if;
F := BaseRing(R); S := PolynomialRing(F, 3); x := S.1; y := S.2; z := S.3;
fhom := S ! (z^(Degree(f)) * Evaluate(f, [ x/z, y/z ]));
return Curve(Scheme(ProjectiveSpace(S), fhom));

end intrinsic;


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


intrinsic ChangeRingCurve(X::Crv, phi::.) -> Crv
{Changes base by a specific homomorphism.}

if Type(X) eq CrvHyp then
    fK, hK := HyperellipticPolynomials(X);
    L := Codomain(phi); S := PolynomialRing(L);
    if fK eq 0 then
        fL := S ! 0;
    else
        fL := &+[ phi(Coefficient(fK, d))*S.1^d : d in [0..Degree(fK)] ];
    end if;
    if hK eq 0 then
        hL := S ! 0;
    else
        hL := &+[ phi(Coefficient(hK, d))*S.1^d : d in [0..Degree(hK)] ];
    end if;
    return HyperellipticCurve(fL, hL);

elif Type(X) eq CrvPln then
    FK := DefiningPolynomials(X)[1];
    L := Codomain(phi); S := PolynomialRing(L, 3);
    FL := &+[ phi(MonomialCoefficient(FK, mon))*Monomial(S, Exponents(mon)) : mon in Monomials(FK) ];
    return PlaneCurve(FL);
end if;

end intrinsic;
