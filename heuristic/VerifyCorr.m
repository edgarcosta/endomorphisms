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
    return P;
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
            rt := Roots(fac, L)[1][1];
            XL := ChangeRing(X, L);
            if z0 ne 0 then
                P := XL ! [ n0*rt + x0, rt + y0, z0 ];
            else
                P := XL ! [ n0*rt + x0, y0, rt + z0 ];
            end if;
            if not IsWeierstrassPlace(Place(P)) then
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
        rt := Roots(h(f), L)[1][1];
        XL := ChangeRing(X, L);
        P := XL ! [ n0, rt, 1 ];
        if not IsWeierstrassPlace(Place(P)) then
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

K := Parent(P[1]); L := Parent(Q[1]);
M, phiK, phiL := RelativeCompositum(K, L);
XM := ChangeRing(X, M); YM := ChangeRing(Y, M);
PM := [ phiK(c) : c in P ]; PM := XM ! P;
QM := [ phiL(c) : c in Q ]; QM := YM ! Q;
AM := ChangeRing(A, M);

if (#Rows(AM) eq #Rows(Transpose(AM))) and IsScalar(AM) then
    return true, "Scalar: OK";
elif Genus(Y) eq 1 then
    return CantorMorphismFromMatrixSplit(XM, PM, YM, QM, AM/2);
else
    return "";
    //return DivisorFromMatrixSplit(XM, PM, YM, QM, AM);
end if;

end intrinsic;
