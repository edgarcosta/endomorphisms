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


function SmallBasePointHyp(X : Bound := 2^10, NW := false)

K := BaseRing(X);
f, h := HyperellipticPolynomials(X);

/* Elliptic case: */
if Genus(X) eq 1 and Degree(f) eq 3 then
    XK := ChangeRing(X, K);
    P := [ K ! 1, K ! 0, K ! 0 ];
    return XK ! P;
end if;

/* Small point over the base field (only QQ permitted for now): */
if Type(BaseRing(X)) eq FldRat then
    g := 4*f + h^2;
    lcm := LCM([ Denominator(c) : c in Coefficients(g) ]);
    a, b := SquareFreeFactorization(lcm);
    g *:= (a*b)^2;
    Y := HyperellipticCurve(g);
    Qs := RationalPoints(Y : Bound := Bound);
    if NW then
        Qs := [ Q : Q in Qs | not IsWeierstrassPlace(Place(Q)) ];
    end if;
    Qs_inf := [ Q : Q in Qs | Q[3] eq 0 ];
    Qs_fin := [ Q : Q in Qs | Q[3] ne 0 ];

    if #Qs_inf ne 0 then
        Q := Qs_inf[1];
        h0 := Coefficient(h, Degree(g) div 2);
        P := [ Q[1], (Q[2] - h0)/(2*a*b), Q[3] ];
        P := [ K ! P[1], K ! P[2], K ! P[3] ];
        XK := ChangeRing(X, K);
        return XK ! P;
    end if;

    if #Qs_fin ne 0 then
        Hts := [ Maximum([ Height(c) : c in Eltseq(Q) ]) : Q in Qs_fin ];
        min, ind := Minimum(Hts);
        Q := Qs_fin[ind];
        h0 := Evaluate(h, Q[1]);
        P := [ Q[1], (Q[2] - h0)/(2*a*b), Q[3] ];
        P := [ K ! P[1], K ! P[2], K ! P[3] ];
        XK := ChangeRing(X, K);
        return XK ! P;
    end if;
end if;

g := 4*f + h^2; Y := HyperellipticCurve(g);
d := Degree(g);

/* Infinite patch: */
if IsEven(d) then
    e := d div 2;
    g0 := Coefficient(g, d);
    h0 := Coefficient(h, e);
    R<t> := PolynomialRing(K);
    /* L is absolute over base field of the curve */
    L, hKL := ExtendNumberFieldExtra(t^2 - g0);
    Q := [ 1, Roots(t^2 - g0, L)[1][1], 0 ];
    P := [ Q[1], (Q[2] - h0)/2, Q[3] ];
    XL := ChangeRingCurve(X, hKL);
    return XL ! P;
end if;

/* Finite patch: */
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
    L, hKL := ExtendNumberFieldExtra(t^2 - g0);
    Q := [ L ! n0, Roots(t^2 - g0, L)[1][1], 1 ];
    h0 := Evaluate(h, Q[1]);
    P := [ Q[1], (Q[2] - h0)/2, Q[3] ];
    XL := ChangeRingCurve(X, hKL);
    return XL ! P;
end if;
error "All cases in SmallBasePointHyp fell through";

end function;


function SmallBasePointPlane(X : Bound := 2^10, NW := false)

K := BaseRing(X);
Ps := RationalPoints(X);
if NW then
    Ps := [ P : P in Ps | not IsWeierstrassPlace(Place(P)) ];
end if;
if #Ps ne 0 then
    Hts := [ Maximum([ Height(c) : c in Eltseq(P) ]) : P in Ps ];
    min, ind := Minimum(Hts);
    P := Ps[ind];
    return X ! Eltseq(P);
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
        for tup in Fac do
            fac := tup[1];
            L, hKL := ExtendNumberFieldExtra(fac);
            rt := Roots(fac, L)[1][1];
            if z0 ne 0 then
                P := [ n0*rt + x0, rt + y0, z0 ];
            else
                P := [ n0*rt + x0, y0, rt + z0 ];
            end if;
            XL := ChangeRingCurve(X, hKL);
            if NW then
                if not IsWeierstrassPlace(Place(XL ! P)) then
                    return P;
                end if;
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
        L, hKL := ImproveFieldExtra(NumberFieldExtra(tup[1]));
        rt := Roots(hRS(f), L)[1][1];
        P := [ n0, rt, 1 ];
        XL := ChangeRingCurve(X, hKL);
        if NW then
            if not IsWeierstrassPlace(Place(XL ! P)) then
                return XL ! P;
            end if;
        end if;
    end for;
    n0 +:= 1;
end while;

end function;


intrinsic SmallBasePoint(X::Crv : Bound := 2^10, NW := NW) -> SeqEnum
{Given a curve X, returns a point on X over a small extension of its base
field, which we can ask to be a non-Weierstrass point.}

if Type(X) eq CrvHyp then
    return SmallBasePointHyp(X : Bound := Bound, NW := NW);
elif Type(X) eq CrvPln then
    return SmallBasePointPlane(X : Bound := Bound, NW := NW);
end if;
error "Not implemented for general curves yet";

end intrinsic;


intrinsic Correspondence(A::., P::SeqEnum, Q::SeqEnum) -> .
{Given curves X and Y with non-Weierstrass points P and Q respectively, finds a
correspondence with tangent representation A if it exists.}

X := Curve(P); Y := Curve(Q);
/* Change everything to common extension: */
// TODO: First take extension for Q, then over that field extension for P, to
// avoid interference, that is, field incompatibilities
KA := Parent(A[1,1]); KP := Parent(P[1]); KQ := Parent(Q[1]);
L, phis := CompositumExtra([* KA, KP, KQ *]);
phiA, phiP, phiQ := Explode(phis);

/* Change to common base */
XL := ChangeRingCurve(X, phiP);
YL := ChangeRingCurve(Y, phiQ);
PL := [ phiP(c) : c in Eltseq(P) ]; PL := XL ! PL;
QL := [ phiQ(c) : c in Eltseq(Q) ]; QL := YL ! QL;
AL := Matrix(L, [ [ phiA(c) : c in Eltseq(row) ] : row in Rows(A) ]);

/* Actual work */
if (#Rows(AL) eq #Rows(Transpose(AL))) and IsScalar(AL) then
    return true, "Scalar: OK for now";
elif Genus(Y) eq 1 then
    test, fs := CantorFromMatrixAmbientSplit(AL, PL, QL);
    if test and (not CorrespondenceVerifyG1(AL, PL, QL)) then
        error "Pullback incorrect";
    end if;
    return test, fs;
else
    return DivisorFromMatrixAmbientSplit(AL, PL, QL);
end if;

end intrinsic;


intrinsic CorrespondenceVerifyG1(A::., fs::SeqEnum, P::Pt, Q::Pt : CheckDegree := false) -> BoolElt
{Returns whether the morphism defined by fs indeed corresponds to the tangent
representation A.}

X := Curve(P); Y := Curve(Q);
F := BaseRing(Parent(fs[1]));
/* Check that the answer is a projection: */
R<x,y> := PolynomialRing(F, 2);
K := FieldOfFractions(R);
fX := R ! DefiningEquation(AffinePatch(X, 1)); fY := R ! DefiningEquation(AffinePatch(Y, 1));
IX := ideal<R | fX>;
if not R ! Numerator(K ! Evaluate(R ! fY, fs)) in IX then
    return false;
end if;

/* Only use if the degree is small, if at all */
if CheckDegree then
    AX := AffinePatch(X, 1); AY := AffinePatch(Y, 1);
    KX := FunctionField(AX); KY := FunctionField(AY);
    m := map<AX -> AY | fs >;
    print "Degree:", Degree(ProjectiveClosure(m));
end if;

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
