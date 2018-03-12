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


intrinsic PlaneCurve(F::RngMPolElt) -> Crv
{Returns the plane curve defined by F.}

R := Parent(F); gens := GeneratorsSequence(R);
if #gens eq 3 then
    return Curve(Scheme(ProjectiveSpace(R), F));
end if;
S := PolynomialRing(BaseRing(R), 3);
x := S.1; y := S.2; z := S.3;
Fhom := S ! (z^(Degree(F)) * Evaluate(F, [ x/z, y/z ]));
return Curve(Scheme(ProjectiveSpace(S), Fhom));

end intrinsic;


intrinsic CurveType(X::Crv) -> MonStgElt
{Returns a string that describes the type of curve that X belongs to, which is
one of "hyperelliptic", "plane" and "general".}

if Type(X) eq CrvHyp then
    return "hyperelliptic";
elif Type(X) eq CrvPln then
    return "plane";
else
    return "general";
end if;

end intrinsic;


intrinsic EmbedCurveEquations(X::Crv, prec::RngIntElt) -> MonStgElt
{Returns the defining equations of X base changed to CC to precision prec,
using the infinite place of the base ring of X.}

if Type(X) eq CrvHyp then
    f, h := HyperellipticPolynomials(X);
    return EmbedAsComplexPolynomials([ f, h ], prec);
elif Type(X) eq CrvPln then
    F := DefiningPolynomial(X);
    return EmbedAsComplexPolynomials([ F ], prec);
end if;
error "Function not available for general curves";

end intrinsic;


intrinsic ChangeRingCurve(X::Crv, phi::.) -> Crv
{Returns X base changed by a morphism phi whose domain is the base ring of X.}
/* ChangeRing needs automatic coercion... */

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
error "Function not available for general curves";

end intrinsic;
