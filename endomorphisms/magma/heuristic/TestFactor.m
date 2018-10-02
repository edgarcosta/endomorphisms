/***
 *  Test if given elliptic curve is a factor of the Jacobian
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic TestEllipticFactor(X::Crv, E::Crv, F::Fld : prec := 300) -> BoolElt
{Given a curve X and an elliptic curve E, determines whether E is a factor of
the Jacobian of X. Returns a map from X to E if this is the case.}

P := PeriodMatrix(X : prec := prec); Q := PeriodMatrix(E : prec := prec);
HomRep := GeometricHomomorphismRepresentation(P, Q, F);
if #HomRep eq 0 then
    return false;
end if;
gen0, d0 := MorphismOfSmallDegree(HomRep, F);
return true, [* gen0, d0 *];

end intrinsic;


intrinsic MorphismOfSmallDegree(HomRep::., F::Fld : Bound := 10) -> List, RngIntElt
{Gives a morphism of small degree from the Jacobian corresponding to P to that
corresponding to Q. The third argument F is the base field used.}

g := #Rows(P);
HomRep := GeometricHomomorphismRepresentation(P, Q, F);
Rs := [ rep[2] : rep in HomRep ];

D := [-Bound..Bound];
M1 := ChangeRing(StandardSymplecticMatrix(1), Rationals());
Mg := ChangeRing(StandardSymplecticMatrix(g), Rationals());
CP := CartesianPower(D, #Rs);
d0 := Infinity();
for tup in CP do
    R := &+[ tup[i] * Rs[i] : i in [1..#Rs] ];
    C := R*Mg*Transpose(R)*M1^(-1);
    test := IsScalar(C);
    d := Abs(C[1,1]);
    if (not IsZero(R)) and d lt d0 then
        d0 := d;
        tup0 := tup;
    end if;
end for;

A0 := &+[ tup0[i]*HomRep[i][1] : i in [1..#HomRep] ];
R0 := &+[ tup0[i]*HomRep[i][2] : i in [1..#HomRep] ];
gen0 := [* A0, R0 *];
return gen0, d0;

end intrinsic;


intrinsic EllipticCMCurve(D::RngIntElt : prec := 1000) -> BoolElt
{Determines principal curve with CM by D.}

QQ := RationalsExtra(prec); CC := QQ`CC; RR := RealField(CC);
CCLarge := ComplexFieldExtra(prec + 100);
if D mod 4 ne 0 then
    tau := (Sqrt(CCLarge ! D) + 1)/2;
else
    tau := Sqrt(CCLarge ! D)/2;
end if;
jCC := CC ! jInvariant(tau);

K, js := NumberFieldExtra(jCC, QQ); j := js[1];
E := EllipticCurveFromjInvariant(j); E := WeierstrassModel(E);

if Type(K) eq FldRat then
    return MinimalModel(E);
end if;
cs := Coefficients(E); a := cs[4]; b := cs[5];
da := Denominator(a); db := Denominator(b);
Fa := Factorization(da); Fb := Factorization(db);
psa := [ tup[1] : tup in Fa ]; psb := [ tup[1] : tup in Fb ];
ps := Set(psa cat psb);
lambda := &*[ p^(Maximum(Ceiling(Valuation(da, p)/2), Ceiling(Valuation(da, p)/3))) : p in ps ];
a *:= lambda^2; b *:= lambda^3;
R<x> := PolynomialRing(K);
f := x^3 + a*x + b;
E := HyperellipticCurve(f);
return E;

end intrinsic;
