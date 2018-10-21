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
{Returns the plane curve defined by F, which can be given affinely or
projectively.}

R := Parent(F); gens := GeneratorsSequence(R);
if #gens eq 3 then
    return Curve(Scheme(ProjectiveSpace(R), F));
end if;
S := PolynomialRing(BaseRing(R), 3);
x := S.1; y := S.2; z := S.3;
Fhom := S ! (z^(Degree(F)) * Evaluate(F, [ x/z, y/z ]));
return Curve(Scheme(ProjectiveSpace(S), Fhom));

end intrinsic;


intrinsic HyperellipticCurveExtra(f::RngUPolElt, h::RngUPolElt, prec::RngIntElt) -> Crv
{Returns the hyperelliptic curve defined by F, which can be given affinely or
projectively.}

QQ := RationalsExtra(prec); RQQ := PolynomialRing(QQ);
return HyperellipticCurve(RQQ ! f, RQQ ! h);

end intrinsic;


intrinsic PlaneCurveExtra(F::RngMPolElt, prec::RngIntElt) -> Crv
{Returns the plane curve defined by F, which can be given affinely or
projectively.}

QQ := RationalsExtra(prec); RQQ := PolynomialRing(QQ, #GeneratorsSequence(Parent(F)));
return PlaneCurve(RQQ ! F);

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


intrinsic EmbedCurveEquations(X::Crv) -> MonStgElt
{Returns the defining equations of X base changed to CC to precision prec,
using the infinite place of the base ring of X.}

if Type(X) eq CrvHyp or Type(X) eq CrvEll then
    f, h := HyperellipticPolynomials(X);
    return EmbedAtInfinitePlace([ f, h ]);
elif Type(X) eq CrvPln then
    F := DefiningPolynomial(X);
    return EmbedAtInfinitePlace([ F ]);
end if;
error "Function not available for general curves";

end intrinsic;


intrinsic ChangeRingCurve(X::Crv, h::.) -> Crv
{Returns X base changed by a morphism h whose domain is the base ring of X.}

if Type(X) eq CrvHyp then
    fK, hK := HyperellipticPolynomials(X);
    L := Codomain(h); S := PolynomialRing(L);
    if fK eq 0 then
        fL := S ! 0;
    else
        fL := &+[ h(Coefficient(fK, d))*S.1^d : d in [0..Degree(fK)] ];
    end if;
    if hK eq 0 then
        hL := S ! 0;
    else
        hL := &+[ h(Coefficient(hK, d))*S.1^d : d in [0..Degree(hK)] ];
    end if;
    return HyperellipticCurve(fL, hL);

elif Type(X) eq CrvPln then
    FK := DefiningPolynomials(X)[1];
    L := Codomain(h); S := PolynomialRing(L, 3);
    FL := &+[ h(MonomialCoefficient(FK, mon))*Monomial(S, Exponents(mon)) : mon in Monomials(FK) ];
    return PlaneCurve(FL);
end if;
error "Function not available for general curves";

end intrinsic;
