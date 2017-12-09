/***
 *  Verifying correspondences
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic NonWeierstrassBasePointHyperelliptic(X::Crv, K::Fld : Bound := 2^10) -> SeqEnum
{Returns a non-Weierstrass point over a small extension of K.}

f, h := HyperellipticPolynomials(X);

/* Elliptic case: */
if Genus(X) eq 1 and Degree(f) eq 3 then
    return [K ! 1, K ! 0, K ! 0];
end if;

/* Small point over the base field (only QQ permitted for now): */
if Type(BaseRing(X)) eq FldRat then
    g := 4*f + h^2;
    lcm := LCM([ Denominator(c) : c in Coefficients(g) ]);
    a, b := SquareFreeFactorization(lcm);
    g *:= (a*b)^2;
    Y := HyperellipticCurve(g);
    Qs := RationalPoints(Y : Bound := Bound);
    Qs_nW := [ Q : Q in Qs | not IsWeierstrassPlace(Place(Q)) ];
    Qs_nW_inf := [ Q : Q in Qs_nW | Q[3] eq 0 ];
    Qs_nW_fin := [ Q : Q in Qs_nW | Q[3] ne 0 ];

    if #Qs_nW_inf ne 0 then
        Q := Qs_nW_inf[1];
        h0 := Coefficient(h, Degree(g) div 2);
        P := [ Q[1], (Q[2] - h0)/(2*a*b), Q[3] ];
        P := [ K ! P[1], K ! P[2], K ! P[3] ];
        return P;
    end if;

    if #Qs_nW_fin ne 0 then
        Hts := [ Maximum([ Height(c) : c in Eltseq(Q) ]) : Q in Qs_nW_fin ];
        min, ind := Minimum(Hts);
        Q := Qs_nW_fin[ind];
        h0 := Evaluate(h, Q[1]);
        P := [ Q[1], (Q[2] - h0)/(2*a*b), Q[3] ];
        P := [ K ! P[1], K ! P[2], K ! P[3] ];
        return P;
    end if;
end if;

g := 4*f + h^2; Y := HyperellipticCurve(g);
d := Degree(g);

/* Non-Weierstrass point in infinite patch: */
if IsEven(d) then
    e := d div 2;
    g0 := Coefficient(g, d);
    h0 := Coefficient(h, e);
    R<t> := PolynomialRing(K);
    L := RelativeSplittingField(t^2 - g0);
    L := ClearFieldDenominator(L);
    Q := [ 1, Roots(t^2 - g0, L)[1][1], 0 ];
    P := [ Q[1], (Q[2] - h0)/2, Q[3] ];
    return P;
end if;

/* Non-Weierstrass point in finite patch: */
if IsOdd(d) then
    n0 := 0;
    while true do
        g0 := Evaluate(g, n0);
        if g0 ne 0 then
            break;
        end if;
        n0 +:= 1;
    end while;
    R<t> := PolynomialRing(K);
    L := RelativeSplittingField(t^2 - g0);
    L := ClearFieldDenominator(L);
    Q := [ L ! n0, Roots(t^2 - g0, L)[1][1], 1 ];
    h0 := Evaluate(h, Q[1]);
    P := [ Q[1], (Q[2] - h0)/2, Q[3] ];
    return P;
end if;

error "All cases in NonWeierstrassBasePointHyperelliptic fell through";

end intrinsic;


intrinsic NonWeierstrassBasePointPlane(X::Crv, K::Fld : Bound := 2^10) -> SeqEnum
{Returns a non-Weierstrass point over a small extension of K.}

Ps := RationalPoints(X);
Ps_nW := [ P : P in Ps | not IsWeierstrassPlace(Place(P)) ];
if #Ps_nW ne 0 then
    Hts := [ Maximum([ Height(c) : c in Eltseq(P) ]) : P in Ps_nW ];
    min, ind := Minimum(Hts);
    L := K;
    P := Ps_nW[ind];
    return Eltseq(P);
end if;

f := DefiningPolynomial(X);
R<x,y,z> := PolynomialRing(K, 3);
S<t> := PolynomialRing(K);

/* If there is a rational point, then take lines through it: */
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
        for tup in Fac do
            fac := tup[1];
            L := NumberField(fac);
            L := ClearFieldDenominator(L);
            rt := Roots(fac, L)[1][1];
            if z0 ne 0 then
                P := [ n0*rt + x0, rt + y0, z0 ];
            else
                P := [ n0*rt + x0, y0, rt + z0 ];
            end if;
            XL := ChangeRing(X, L);
            if not IsWeierstrassPlace(Place(XL ! P)) then
                return P;
            end if;
        end for;
        n0 +:= 1;
    end while;
end if;

/* Otherwise lines through (0,0,1): */
n0 := 0;
while true do
    h := hom< R -> S | [ n0, t, 1 ]>;
    Fac := Factorization(h(f));
    for tup in Fac do
        L := NumberField(tup[1]);
        L := ClearFieldDenominator(L);
        rt := Roots(h(f), L)[1][1];
        P := [ n0, rt, 1 ];
        XL := ChangeRing(X, L);
        if not IsWeierstrassPlace(Place(XL ! P)) then
            return P;
        end if;
    end for;
    n0 +:= 1;
end while;

end intrinsic;


intrinsic NonWeierstrassBasePoint(X::Crv, K::Fld : Bound := 2^10) -> SeqEnum
{Returns a non-Weierstrass point over a small extension of K.}

if Type(X) eq CrvHyp then
    return NonWeierstrassBasePointHyperelliptic(X, K : Bound := Bound);
elif Type(X) eq CrvPln then
    return NonWeierstrassBasePointPlane(X, K : Bound := Bound);
else
    error "Not implemented for general curves yet";
end if;

end intrinsic;


intrinsic Correspondence(X::Crv, P::SeqEnum, Y::Crv, Q::SeqEnum, A::.) -> .
{Gives certificate.}

/* Change everything to common extension: */
KY := BaseRing(Y); KP := Parent(P[1]); KQ := Parent(Q[1]); KA := Parent(A[1,1]);
L, phis := RelativeCompositum([* KY, KP, KQ, KA *]);
phiY, phiP, phiQ, phiA := Explode(phis);
XL := ChangeRing(X, L);
YL := ChangeRingCurve(Y, phiY);
PL := [ phiP(c) : c in Eltseq(P) ]; PL := XL ! PL;
QL := [ phiQ(c) : c in Eltseq(Q) ]; QL := YL ! QL;
AL := Matrix(L, [ [ phiA(c) : c in Eltseq(row) ] : row in Rows(A) ]);

/* Make absolute: */
M := AbsoluteField(L);
XM := ChangeRing(XL, M);
YM := ChangeRing(YL, M);
PM := [ M ! c : c in Eltseq(PL) ]; PM := XM ! PM;
QM := [ M ! c : c in Eltseq(QL) ]; QM := YM ! QM;
AM := Matrix(M, [ [ M ! c : c in Eltseq(row) ] : row in Rows(AL) ]);

if (#Rows(AM) eq #Rows(Transpose(AM))) and IsScalar(AM) then
    return true, "Scalar: OK";
elif Genus(Y) eq 1 then
    test, fs := CantorFromMatrixSplit(XM, PM, YM, QM, AM : UpperBound := 50);
    if test and (not CorrespondenceVerifyG1(XM, PM, YM, QM, AM, fs)) then
        error "Pullback incorrect";
    end if;
    return test, fs;
else
    return DivisorFromMatrixSplit(XM, PM, YM, QM, AM : UpperBound := 50);
end if;

end intrinsic;


intrinsic CorrespondenceVerifyG1(X::Crv, P::Pt, Y::Crv, Q::Pt, A::., fs::SeqEnum) -> .
{Checks alleged properties of fs.}

F := BaseRing(Parent(fs[1]));
/* Check that the answer is a projection: */
R<x,y> := PolynomialRing(F, 2);
K := FieldOfFractions(R);
fX := R ! DefiningEquation(AffinePatch(X, 1)); fY := R ! DefiningEquation(AffinePatch(Y, 1));
IX := ideal<R | fX>;
if not R ! Numerator(K ! Evaluate(R ! fY, fs)) in IX then
    return false;
end if;

/*
AX := AffinePatch(X, 1); AY := AffinePatch(Y, 1);
KX := FunctionField(AX); KY := FunctionField(AY);
m := map<AX -> AY | fs >;
print "Degree:", Degree(ProjectiveClosure(m));
*/

/* Check that the action on differentials is correct: */
fX := R ! DefiningEquation(AffinePatch(X, 1));
fY := R ! DefiningEquation(AffinePatch(Y, 1));
dx := K ! 1;
dy := K ! -Derivative(fX, 1) / Derivative(fX, 2);
ev := ((K ! Derivative(fs[1], 1))*dx + (K ! Derivative(fs[1], 2))*dy) / (K ! (2*fs[2]));
if X`is_hyperelliptic then
    if not R ! Numerator(K ! (ev - &+[ A[1,i]*x^(i - 1) : i in [1..Genus(X)] ]/(Derivative(fX, 2)/MonomialCoefficient(fX, y^2)))) in IX then
        return false;
    end if;
elif X`is_planar then
    mults := [ x, y, 1 ];
    print "Correct pullback?", R ! Numerator(K ! (ev - &+[ A[1,i]*mults[i] : i in [1..Genus(X)] ]/Derivative(fX, 2))) in IX;
end if;
return true;

end intrinsic;
