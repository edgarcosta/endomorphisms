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


intrinsic TestEllipticFactor(X::Crv, E::Crv) -> BoolElt
{Given a curve X and an elliptic curve E, determines whether E is a factor of
the Jacobian of X. Returns a map from X to E if this is the case.}

return 0;

end intrinsic;


intrinsic EllipticCMCurve(D::RngIntElt : prec := 1000) -> BoolElt
{Determines principal curve with CM by D.}

CC := ComplexFieldExtra(prec); RR := RealField(CC);
CCLarge := ComplexFieldExtra(prec + 100);
if D mod 4 ne 0 then
    tau := (Sqrt(CCLarge ! D) + 1)/2;
else
    tau := Sqrt(CCLarge ! D)/2;
end if;
jCC := CC ! jInvariant(tau);

QQ := RationalsExtra();
p := RelativeMinimalPolynomial(jCC, QQ);
R<t> := PolynomialRing(QQ);
if Degree(p) eq 1 then
    K := RationalsExtra();
    j := Roots(p, K)[1][1];
else
    K<r> := NumberFieldExtra(p);
    j := K.1;
end if;

E := EllipticCurveFromjInvariant(j); E := WeierstrassModel(E);
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


intrinsic MorphismOfSmallDegree(P::., Q::., F::. : Bound := 10) -> .
{Wat it sez on the tin}

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
ACC0 := &+[ tup0[i]*HomRep[i][3] : i in [1..#HomRep] ];
gen0 := [* A0, R0, ACC0 *];
return gen0, d0;

end intrinsic;


intrinsic MorphismOfSmallDegreePartial(P1::., P2::. : Bound := 10) -> .
{Wat it sez on the tin}

g1 := #Rows(P1); g2 := #Rows(P2);
HomRep := GeometricHomomorphismRepresentationCC(P1, P2);
Rs := [ rep[2] : rep in HomRep ];

D := [-Bound..Bound];
M1 := ChangeRing(StandardSymplecticMatrix(g1), Rationals());
M2 := ChangeRing(StandardSymplecticMatrix(g2), Rationals());
CP := CartesianPower(D, #Rs);
d0 := Infinity();
for tup in CP do
    R := &+[ tup[i] * Rs[i] : i in [1..#Rs] ];
    C := R*M1*Transpose(R)*M2^(-1);
    test := IsScalar(C);
    if test then
        d := Abs(C[1,1]);
        if (not IsZero(R)) and d lt d0 then
            d0 := d;
            tup0 := tup;
        end if;
    end if;
end for;

if d0 eq Infinity() then
    return 0, d0;
end if;
ACC0 := &+[ tup0[i]*HomRep[i][1] : i in [1..#HomRep] ];
R0 := &+[ tup0[i]*HomRep[i][2] : i in [1..#HomRep] ];
gen0 := [* ACC0, R0 *];
return gen0, d0;

end intrinsic;
