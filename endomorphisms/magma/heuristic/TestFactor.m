/***
 *  Test if given (elliptic?) curve is a factor of the Jacobian
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic EllipticCMCurve(D::RngIntElt, F::Fld) -> BoolElt
{Determines principal curve with CM by D over some extension of the base F.}

CC := F`CC; RR := RealField(CC);
CCLarge := ComplexFieldExtra(Precision(CC) + 100);
if -D mod 4 eq 1 then
    tau := (Sqrt(CCLarge ! -D) + 1)/2;
elif -D mod 4 eq 0 then
    tau := Sqrt(CCLarge ! -D)/2;
else
    error "D is not the discriminant of a quadratic order";
end if;
jCC := CC ! jInvariant(tau);

K, js := NumberFieldExtra([ jCC ], F); j := js[1];
E := EllipticCurveFromjInvariant(j); E := WeierstrassModel(E);

if Type(K) eq FldRat then
    E0 := MinimalModel(E);
    return E0;
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


function MorphismOfSmallDegreePQ(P, Q, K : Bound := 3)
// Gives a morphism of small degree from the Jacobian corresponding to P to
// that corresponding to Q. The third argument K is the base field used.

gX := #Rows(P); gY := #Rows(Q);
HomRep, hKM := GeometricHomomorphismRepresentation(P, Q, K);
Rs := [ rep[2] : rep in HomRep ];

D := [-Bound..Bound];
MgX := StandardSymplecticMatrix(gX);
MgY := StandardSymplecticMatrix(gY);

CP := CartesianPower(D, #Rs); d0 := Infinity();
for tup in CP do
    R := &+[ tup[i] * Rs[i] : i in [1..#Rs] ];
    C := R*MgX*Transpose(R)*MgY^(-1);
    /* TODO: This is probably only relevant if we want a true morphism.
     *       See what the degree of the correspondence is otherwise. */
    if IsScalar(C) then
        d := Abs(C[1,1]);
        if (not IsZero(R)) and d lt d0 then
            d0 := d;
            tup0 := tup;
        end if;
    end if;
end for;
if d0 eq Infinity() then
    error "No morphism found";
end if;

A0 := &+[ tup0[i]*HomRep[i][1] : i in [1..#HomRep] ];
M := Codomain(hKM);
L, hLM, seq := SubfieldExtra(M, Eltseq(A0) cat [ hKM(K.1) ] );
if IsQQ(K) then
    hKL := hom< K -> L | >;
else
    hKL := hom< K -> L | seq[#seq] >;
end if;

A0 := CoerceToSubfieldMatrix(A0, M, L, hLM);
R0 := &+[ tup0[i]*HomRep[i][2] : i in [1..#HomRep] ];
return [* A0, R0 *], d0, hKL;

end function;


intrinsic MorphismOfSmallDegree(X::Crv, Y::Crv : Bound := 3) -> .
{Gives a morphism of small degree from X to Y and the base extension of Y to the required field. We allow Y to be defined over an extension of the base field of X.}

F := BaseRing(X); K := BaseRing(Y);
P := PeriodMatrix(X); Q := PeriodMatrix(Y);
mor, d0, hKL := MorphismOfSmallDegreePQ(P, Q, K : Bound := Bound);
return [* ChangeRingCurve(Y, hKL), mor *], d0;

end intrinsic;
