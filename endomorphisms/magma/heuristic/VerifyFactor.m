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


intrinsic MorphismOfSmallDegreeHeuristic(P::., Q::., F::Fld : Bound := 10) -> List, RngIntElt
{Gives a morphism of small degree from the Jacobian corresponding to P to that
corresponding to Q. The third argument F is the base field used.}

gX := #Rows(P); gY := #Rows(Q);
HomRep := GeometricHomomorphismRepresentation(P, Q, F);
Rs := [ rep[2] : rep in HomRep ];

D := [-Bound..Bound];
MgX := ChangeRing(StandardSymplecticMatrix(gX), Rationals());
MgY := ChangeRing(StandardSymplecticMatrix(gY), Rationals());
CP := CartesianPower(D, #Rs);
d0 := Infinity();
for tup in CP do
    R := &+[ tup[i] * Rs[i] : i in [1..#Rs] ];
    C := R*MgX*Transpose(R)*MgY^(-1);
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


intrinsic MorphismOfSmallDegreeHeuristic(X::Crv, Y::Crv : Bound := 10) -> .
{Gives a morphism of small degree from X to Y.}

F := BaseRing(X);
P := PeriodMatrix(X); Q := PeriodMatrix(Y);
return MorphismOfSmallDegree(P, Q, F);

end intrinsic;


intrinsic EllipticCMCurve(D::RngIntElt : prec := 10000) -> BoolElt
{Determines principal curve with CM by D.}

QQ := RationalsExtra(prec); CC := QQ`CC; RR := RealField(CC);
CCLarge := ComplexFieldExtra(prec + 100);
if -D mod 4 ne 0 then
    tau := (Sqrt(CCLarge ! -D) + 1)/2;
else
    tau := Sqrt(CCLarge ! -D)/2;
end if;
jCC := CC ! jInvariant(tau);

K, js := NumberFieldExtra([ jCC ], QQ); j := js[1];
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
