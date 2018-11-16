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

declare attributes Crv : plane_model, period_matrix, geo_endo_rep, base_endo_rep;


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


function EmbedCurveEquations(X)
// Returns the defining equations of X base changed to CC to precision prec,
// using the infinite place of the base ring of X.

if Type(X) eq CrvHyp or Type(X) eq CrvEll then
    f, h := HyperellipticPolynomials(X);
    return EmbedAtInfinitePlacePolynomials([ f, h ]);
elif Type(X) eq CrvPln then
    F := DefiningPolynomial(X);
    return EmbedAtInfinitePlacePolynomials([ F ]);
end if;
error "Function not available for general curves";

end function;


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


function CurveDescriptionHyperelliptic(X, F)

if Genus(X) eq 1 then
    desc := "ell";

    K := BaseRing(X); F := BaseRing(K);
    K_seq := FieldDescriptionExtra(K, F);
    field := K_seq;

    f, h := HyperellipticPolynomials(X);
    f_seq := Eltseq(f); h_seq := Eltseq(h);
    f_seq_seq := [ ElementDescriptionExtra(coeff, F) : coeff in f_seq ];
    h_seq_seq := [ ElementDescriptionExtra(coeff, F) : coeff in h_seq ];
    coeffs := [ f_seq_seq, h_seq_seq ];

else
    desc := "hyp";
    field := [ ];
    coeffs := [ ];
end if;

return [* desc, field, coeffs *];

end function;


function CurveDescriptionPlane(X, F)

desc := "pln";

K := BaseRing(X); F := BaseRing(K);

K_seq := FieldDescriptionExtra(K, F);
field := K_seq;

f := DefiningPolynomials(X)[1];
mons := Monomials(f);
coeffsexps := [ ];
for mon in mons do
    coeff := MonomialCoefficient(f, mon);
    coeff_seq := ElementDescriptionExtra(coeff, F);
    exp := Exponents(mon);
    Append(~coeffsexps, [* coeff_seq, exp *]);
end for;

return [* desc, field, coeffsexps *];

end function;


intrinsic CurveDescription(X::Crv, F::Fld) -> List
{Returns a string description of the curve X over the field F.}

if Type(X) eq CrvHyp then
    return CurveDescriptionHyperelliptic(X, F);
elif Type(X) eq CrvPln then
    return CurveDescriptionPlane(X, F);
else
    error "No description for general curves yet";
end if;

end intrinsic;


/* These functions are only used in Sage */
intrinsic HyperellipticCurveExtra(f::RngUPolElt, h::., prec::RngIntElt) -> Crv
{Returns the hyperelliptic curve over the rationals with precision prec defined
by f and h. Only relevant in the Sage interface.}

QQ := RationalsExtra(prec); RQQ := PolynomialRing(QQ);
return HyperellipticCurve(RQQ ! f, RQQ ! h);

end intrinsic;


intrinsic PlaneCurveExtra(F::RngMPolElt, prec::RngIntElt) -> Crv
{Returns the plane curve over the rationals with precision prec defined by F,
which can be given affinely or projectively. Only relevant in the Sage interface.}

QQ := RationalsExtra(prec); RQQ := PolynomialRing(QQ, #GeneratorsSequence(Parent(F)));
return PlaneCurve(RQQ ! F);

end intrinsic;


intrinsic PlaneModel(X::Crv) -> .
{Returns a plane model for more complicated curves.}

if assigned X`plane_model then
    return X`plane_model;
end if;

if Type(X) eq CrvHyp or Type(X) eq CrvPln then
    X`plane_model := X;
    return X`plane_model;
end if;

F := DefiningPolynomial(AlgorithmicFunctionField(FunctionField(X)));
L := Parent(F); K := BaseRing(L);

R := PolynomialRing(BaseRing(K));
S := PolynomialRing(BaseRing(K), 2);
h := hom< R -> S | S.1 >;
F0 := S ! 0;

coeffs := Coefficients(F);
for i in [1..#coeffs] do
    F0 +:= h(R ! coeffs[i])*S.2^(i - 1);
end for;
X`plane_model := PlaneCurve(F0);
return X`plane_model;

end intrinsic;
