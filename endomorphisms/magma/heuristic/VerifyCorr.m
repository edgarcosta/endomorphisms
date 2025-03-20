/***
 *  Verifying correspondences
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@unipi.it)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

forward SmallBasePointHyp;
forward SmallBasePointPlane;
forward CorrespondenceVerifyG1;


function SmallBasePointHyp(X : Bound := 2^10, NW := false)

K := BaseRing(X);
f, h := HyperellipticPolynomials(X);

/* Elliptic case: */
if Genus(X) eq 1 and Degree(f) eq 3 then
    P := [ 1, 0, 0 ];
    return X ! P, CanonicalInclusionMap(K, K);
end if;

/* Small point over the base field (only QQ permitted for now): */
if Type(K) eq FldRat then
    g := 4*f + h^2;
    lcm := LCM([ Denominator(c) : c in Coefficients(g) ]);
    a, b := SquareFreeFactorization(lcm);
    g *:= (a*b)^2;
    Y := HyperellipticCurve(g);
    Qs := RationalPoints(Y : Bound := Bound);

    /* Exclude Weierstrass points among small points */
    if NW then
        Qs := [ Q : Q in Qs | not IsWeierstrassPlace(Place(Q)) ];
    end if;

    /* Partition any points returned into infinite and finite points */
    Qs_inf := [ Q : Q in Qs | Q[3] eq 0 ];
    Qs_fin := [ Q : Q in Qs | Q[3] ne 0 ];

    /* First try to transform back infinite place */
    if #Qs_inf ne 0 then
        Q := Qs_inf[1];
        h0 := Coefficient(h, Degree(g) div 2);
        P := [ Q[1], (Q[2] - h0)/(2*a*b), Q[3] ];
        P := [ P[1], P[2], P[3] ];
        return X ! P, CanonicalInclusionMap(K, K);
    end if;

    /* Otherwise transform back finite place */
    if #Qs_fin ne 0 then
        Hts := [ Maximum([ Height(c) : c in Eltseq(Q) ]) : Q in Qs_fin ];
        min, ind := Minimum(Hts);
        Q := Qs_fin[ind];
        h0 := Evaluate(h, Q[1]);
        P := [ Q[1], (Q[2] - h0)/(2*a*b), Q[3] ];
        P := [ P[1], P[2], P[3] ];
        return X ! P, CanonicalInclusionMap(K, K);
    end if;
end if;

/* Preparing to deal with general base field */
g := 4*f + h^2; Y := HyperellipticCurve(g);
d := Degree(g);

/* Prefer infinite places: */
if IsEven(d) or not NW then
    e := d div 2;
    g0 := Coefficient(g, 2*e); h0 := Coefficient(h, e);
    R<t> := PolynomialRing(K);

    rts := Roots(t^2 - g0);
    L := K; hKL := CanonicalInclusionMap(K, K);
    if #rts eq 0 then
        L, r, hKL := NumberFieldExtra(t^2 - g0);
        X := ChangeRingCurve(X, hKL);
        rts := [ [ r, 1 ] ];
    end if;
    Q := [ 1, rts[1][1], 0 ];
    P := [ Q[1], (Q[2] - hKL(h0))/2, Q[3] ];
    return X ! P, hKL;
end if;

/* Finite patch: */
if IsOdd(d) then
    n0 := 0;
    /* Next if NW is superfluous */
    if NW then
        while true do
            g0 := Evaluate(g, n0);
            if g0 ne 0 then
                break;
            end if;
            n0 +:= 1;
        end while;
    end if;
    R<t> := PolynomialRing(K);

    rts := Roots(t^2 - g0);
    L := K; hKL := CanonicalInclusionMap(K, K);
    if #rts eq 0 then
        L, r, hKL := NumberFieldExtra(t^2 - g0);
        X := ChangeRingCurve(X, hKL);
        rts := [ [ r, 1 ] ];
    end if;
    Q := [ n0, rts[1][1], 1 ];
    h0 := Evaluate(h, Q[1]);
    P := [ Q[1], (Q[2] - hKL(h0))/2, Q[3] ];
    return X ! P, hKL;
end if;
error "All cases in SmallBasePointHyp fell through";

end function;


function SmallBasePointPlane(X : Bound := 2^10, NW := false)

/* Choose rational point of small height if possible */
K := BaseRing(X);
Ps := RationalPoints(X : Bound := Bound);
if NW then
    Ps := [ P : P in Ps | not IsWeierstrassPlace(Place(P)) ];
end if;
if #Ps ne 0 then
    Hts := [ Maximum([ Height(c) : c in Eltseq(P) ]) : P in Ps ];
    min, ind := Minimum(Hts);
    P := Ps[ind];
    return X ! Eltseq(P), CanonicalInclusionMap(K, K);
end if;

f := DefiningPolynomial(X);
R<x,y,z> := PolynomialRing(K, 3);
S<t> := PolynomialRing(K);

/* If there is a rational point, then take lines through it */
if #Ps ne 0 then
    P0 := Ps[1];
    x0, y0, z0 := Explode(Eltseq(P0));
    n0 := 0;
    while true do
        if z0 ne 0 then
            h := hom< R -> S | [ n0*t + x0, t + y0, z0 ]>;
        else
            h := hom< R -> S | [ n0*t + x0, y0, t + z0 ]>;
        end if;
        Fac := Factorization(h(f));

        /* Factorize and take corresponding extension */
        for tup in Fac do
            fac := tup[1];
            L, rt, hKL := NumberFieldExtra(fac);
            if z0 ne 0 then
                P := [ n0*rt + x0, rt + y0, z0 ];
            else
                P := [ n0*rt + x0, y0, rt + z0 ];
            end if;
            XL := ChangeRingCurve(X, hKL);

            if NW then
                if not IsWeierstrassPlace(Place(XL ! P)) then
                    return XL ! P, hKL;
                end if;
            else
                return XL ! P, hKL;
            end if;
        end for;

        n0 +:= 1;
    end while;
end if;

/* Otherwise lines through (0,0,1): */
n0 := 0;
while true do
    hRS := hom< R -> S | [ n0, t, 1 ]>;
    Fac := Factorization(hRS(f));

    for tup in Fac do
        L, rt, hKL := NumberFieldExtra(tup[1]);
        P := [ n0, rt, 1 ];
        XL := ChangeRingCurve(X, hKL);

        if NW then
            if not IsWeierstrassPlace(Place(XL ! P)) then
                return XL ! P, hKL;
            end if;
        else
            return XL ! P, hKL;
        end if;
    end for;
    n0 +:= 1;
end while;

end function;


intrinsic SmallBasePoint(X::Crv : Bound := 2^10, NW := false) -> SeqEnum
{Given a curve X, returns a point on X over a small extension of its base
field, which we can ask to be a non-Weierstrass point.}

if assigned X`base_point then
    return X`base_point;
end if;

if Type(X) eq CrvHyp then
    P, hFK := SmallBasePointHyp(X : Bound := Bound, NW := NW);
elif Type(X) eq CrvPln then
    P, hFK := SmallBasePointPlane(X : Bound := Bound, NW := NW);
else
    error "Not implemented for general curves yet";
end if;

X`base_point := [* P, hFK *];
return X`base_point;

end intrinsic;


intrinsic Correspondence(X::Crv, Y::Crv, mor::. : P := 0, Q := 0, CheckDegree := false, Al := "Cantor", B := 300) -> .
{Given curves X and Y, finds a correspondence with tangent representation A if it exists.
Base points P and Q can be specified: otherwise these are found automatically over some extension.}

assert Al in [ "Cantor", "Divisor" ];
A := mor[1];
R := mor[2];
KX := BaseRing(X);
KY := BaseRing(Y);
FX, hFKX := InclusionOfBaseExtra(KX);
FY, hFKY := InclusionOfBaseExtra(KY);
require FY eq FX: "different bases for KX and KY";
KA := BaseRing(Parent(A));
FA, hFKA := InclusionOfBaseExtra(KA);
require FA eq FX: "different bases for KX and KA";

/* Use or find point on X */
if Type(P) ne RngIntElt then
    assert not IsWeierstrassPlace(Place(X ! P));
    hKXKP := CanonicalInclusionMap(KX, KX);
    KP := KX;
else
P, hKXKP := Explode(SmallBasePoint(X : NW:=true));
    KP := BaseRing(Curve(P));
end if;

/* Use or find point on Y */
if Type(Q) ne RngIntElt then
    assert not IsWeierstrassPlace(Place(Y ! Q));
    hKYKQ := CanonicalInclusionMap(KY, KY);
    KQ := KY;
else
Q, hKYKQ := Explode(SmallBasePoint(Y : NW:=true));
    KQ := BaseRing(Curve(Q));
end if;
hFKQ := hFKY*hKYKQ;

/* Change to common base field */
KAP, hKAKAP, hKPKAP := CompositumExtra(KA, KP);
KAPQ, hKAPKAPQ, hKQKAPQ := CompositumExtra(KAP, KQ);
hFKAPQ := hFKQ*hKQKAPQ;
hKPKAPQ := hKPKAP*hKAPKAPQ;
hKXKAPQ := hKXKP*hKPKAPQ;
hKYKAPQ := hKYKQ*hKQKAPQ;
hKAKAPQ := hKAKAP*hKAPKAPQ;

XKAPQ := ChangeRingCurve(X, hKXKAPQ);
PKAPQ := XKAPQ ! [ hKPKAPQ(c) : c in Eltseq(P) ];
YKAPQ := ChangeRingCurve(Y, hKYKAPQ);
QKAPQ := YKAPQ ! [ hKQKAPQ(c) : c in Eltseq(Q) ];
AKAPQ := ConjugateMatrix(hKAKAPQ, A);

/* Actual work */
X := XKAPQ; P := PKAPQ; Y := YKAPQ; Q := QKAPQ; A := AKAPQ;
vprint EndoCheck : "";
vprint EndoCheck : "Chosen base point on X:";
vprint EndoCheck : P;
vprint EndoCheck : "Chosen base point on Y:";
vprint EndoCheck : Q;

if (ISA(Type(R), Mtrx) and #Rows(R) eq #Rows(Transpose(R))) and IsScalar(R) then
    return true, "Multiplication by an integer";
else
    if Al eq "Cantor" then
	test, fs := CantorFromMatrixAmbientSplit(X, P, Y, Q, A : B := B );
        if Genus(Y) eq 1 then
            if test and (not CorrespondenceVerifyG1(X, Y, A, fs : CheckDegree := CheckDegree)) then
                error "Pullback incorrect";
            end if;
        end if;
        return test, fs;
    elif Al eq "Divisor" then
        test, D := DivisorFromMatrixAmbientSplit(X, P, Y, Q, A : B := B);
        return test, D;
    end if;
end if;

end intrinsic;


intrinsic CorrespondenceVerifyG1(X::Crv, Y::Crv, A::., fs::SeqEnum : CheckDegree := false) -> BoolElt
{Returns whether the morphism defined by fs indeed corresponds to the tangent
representation A.}

gY := Y`g;
if gY eq 1 then
    return CorrespondenceVerifyG1(X, Y, A, fs : CheckDegree := CheckDegree);
else
    error "No verification algorithm implemented yet in this case";
end if;

end intrinsic;


function CorrespondenceVerifyG1(X, Y, A, fs : CheckDegree := false)
// Returns whether the morphism defined by fs indeed corresponds to the tangent
// representation A.

/* Check that the answer is well-defined */
F := BaseRing(Parent(fs[1]));
R<x,y> := PolynomialRing(F, 2);
K := FieldOfFractions(R);
fX := R ! DefiningEquation(AffinePatch(X, 1));
fY := R ! DefiningEquation(AffinePatch(Y, 1));
IX := ideal<R | fX>;
if not R ! Numerator(K ! Evaluate(R ! fY, fs)) in IX then
    return false;
end if;

/* Only use if the degree is small, if at all */
if CheckDegree then
    AX := AffinePatch(X, 1); AY := AffinePatch(Y, 1);
    KX := FunctionField(AX); KY := FunctionField(AY);
    m := map<AX -> AY | fs >;
    vprint EndoCheck : "";
    vprint EndoCheck : "Degree:";
    vprint EndoCheck : Degree(ProjectiveClosure(m));
end if;

/* Check that the action on differentials is correct: */
fX := R ! DefiningEquation(AffinePatch(X, 1));
fY, hY := HyperellipticPolynomials(Y);
dx := K ! 1;
dy := K ! -Derivative(fX, 1) / Derivative(fX, 2);
ev := ((K ! Derivative(fs[1], 1))*dx + (K ! Derivative(fs[1], 2))*dy) / (K ! (2*fs[2] + Evaluate(hY, fs[1])));
if X`is_hyperelliptic then
    if not R ! Numerator(K ! (ev - &+[ A[1,i]*x^(i - 1) : i in [1..Genus(X)] ]/(Derivative(fX, 2)/MonomialCoefficient(fX, y^2)))) in IX then
        return false;
    end if;
elif X`is_planar then
    mults := [ x, y, 1 ];
    vprint EndoCheck : "";
    vprint EndoCheck : "Correct pullback?";
    vprint EndoCheck : R ! Numerator(K ! (ev - &+[ A[1,i]*mults[i] : i in [1..Genus(X)] ]/Derivative(fX, 2))) in IX;
end if;
return true;

end function;
