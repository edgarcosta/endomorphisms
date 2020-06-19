/***
 *  Determining period matrices
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@unipi.it)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

import "Curves.m": EmbedCurveEquations;
forward PeriodMatrixGH;

/* TODO: Add integration over CC via cheat using QQ (i) */


function RiemannSurfaceMN(gCC, e : Precision := 30)
return SE_Curve(gCC, e : Prec := Precision);
end function;


intrinsic PeriodMatrix(X::Crv) -> ModMatFldElt
{Returns the period matrix of X.}

F := BaseRing(X); CC := Parent(F`iota);
vprint EndoFind : "";
vprint EndoFind : "Calculating period matrix...";
if assigned X`period_matrix then
    vprint EndoFind : "using stored period matrix.";
    return X`period_matrix;
end if;

if CurveType(X) eq "genhyp" then
    q, f := Explode(X`ghpols);
    X`period_matrix := PeriodMatrixGH(q, f);
    return X`period_matrix;
end if;

if CurveType(X) eq "hyp" then
    f, h := HyperellipticPolynomials(X);
    g := (4*f + h^2) / 4;
    gCC := EmbedPolynomialExtra(g);
    RS := RiemannSurfaceMN(gCC, 2 : Precision := Precision(CC));
    P := ChangeRing(RS`BigPeriodMatrix, CC) / 2;
end if;

if CurveType(X) in [ "plane", "gen" ] then
    Y := PlaneModel(X);

    /* Determine correct embedding */
    if not Type(F) eq FldRat then
        sigmas := InfinitePlaces(F);
        for sigmatry in sigmas do
            embF := EmbedExtra(F.1);
            embsigma := CC ! Evaluate(F.1, sigmatry : Precision := Precision(CC));
            embsigmacc := ComplexConjugate(embsigma);
            if Abs(embF - embsigma) lt CC`epscomp then
                found := true; cc := false; sigma := sigmatry; break;
            elif Abs(embF - embsigmacc) lt CC`epscomp then
                found := true; cc := true; sigma := sigmatry; break;
            end if;
        end for;
        assert found;
    end if;

    f := DefiningPolynomial(AffinePatch(Y, 1));
    if Type(F) eq FldRat then
        RS := RiemannSurface(f : Precision := Precision(CC));
        P := ChangeRing(BigPeriodMatrix(RS), CC);
    else
        RS := RiemannSurface(f, sigma : Precision := Precision(CC));
        P := ChangeRing(BigPeriodMatrix(RS), CC);
        /* Conjugate if needed */
        if cc then
            P := Matrix(CC, [ [ ComplexConjugate(c) : c in Eltseq(row) ] : row in Rows(P) ]);
        end if;
    end if;
end if;

/* Storing info and ensuring compatibility with Birkenhake--Lange */
X`riesrf := RS;
g := #Rows(P);
P1 := Submatrix(P, 1,1,   g,g);
P2 := Submatrix(P, 1,g+1, g,g);
P := HorizontalJoin(P2, P1);

/* Take canonical form in genus 3 */
if CurveType(X) eq "plane" and g eq 3 then
    rowsP := Rows(P);
    P := Matrix([ rowsP[3], rowsP[2], rowsP[1] ]);
end if;

X`period_matrix := P;
vprint EndoFind : "done calculating period matrix.";
return X`period_matrix;

end intrinsic;


function PeriodMatrixGH(q, f)
/* Compability to calculate with geometrically hyperelliptic curves */

S := Parent(q); F := BaseRing(S); K := F;
S3 := S; K3 := FieldOfFractions(S3); PP2 := ProjectiveSpace(S3);
S2 := PolynomialRing(F, 2); K2 := FieldOfFractions(S2); PP1 := ProjectiveSpace(S2);
R := PolynomialRing(F); h21 := hom< S2 -> R | [ R.1, 1 ] >;

/* Find point on conic and parametrize */
Q := Conic(PP2, q);
test, P := HasRationalPoint(Q);
if not test then
    /* Find suitable base extension */
    p := Evaluate(q, [ R.1, 0, 1 ]); K := NumberFieldExtra(p);
    S := PolynomialRing(K, 3); h := hom< Parent(q) -> S | [ S.1, S.2, S.3 ] >;
    q := h(q); f := h(f);

    /* Redefine */
    S3 := S; KS3 := FieldOfFractions(S3); PP2 := ProjectiveSpace(S3);
    S2 := PolynomialRing(K, 2); KS2 := FieldOfFractions(S2); PP1 := ProjectiveSpace(S2);
    R := PolynomialRing(K); h21 := hom< S2 -> R | [ R.1, 1 ] >;

    /* Now we have a point! */
    Q := Conic(PP2, q);
    p := Evaluate(q, [ R.1, 0, 1 ]); r := Roots(p)[1][1];
    P := Q ! [ r, 0, 1 ];
end if;
phi := Parametrization(Q, P);
DE := DefiningEquations(phi);
h := hom< Parent(DE[1]) -> S2 | [ S2.1, S2.2 ] >;
DE := [ h(c) : c in DE ];

/* Create hyperelliptic curve */
// TODO: This could give incompatibilities! Not yet though.
// [Update: I no longer understand this comment, but I keep it in just in case.]
g := Evaluate(f, DE);
PP2W := ProjectiveSpace(K, [Degree(f),1,1]);
S3W := CoordinateRing(PP2W);
h23 := hom< S2 -> S3W | [ S3.2, S3.3 ] >;
Y := Curve(PP2W, [ S3W.1^2 - h23(g) ]);

/* Create original curve */
PP3 := ProjectiveSpace(K, [2,1,1,1]);
S4 := CoordinateRing(PP3);
h34 := hom< S3 -> S4 | [ S4.2, S4.3, S4.4 ] >;
eqs := [ h34(q), S4.1^2 - h34(f) ];
X := Curve(PP3, eqs);

/* Isomorphism between the and corresponding pullback */
m := map< Y -> X | [ S3W.1, h23(DE[1]), h23(DE[2]), h23(DE[3]) ] >;
BX := BasisOfHolomorphicDifferentials(X);
BXPB := [ Pullback(m, omega) : omega in BX ];
BY := BasisOfHolomorphicDifferentials(Y);

/* Read off matrix on differentials */
genus := Degree(f) - 1;
rows := [ ];
for omega in BXPB do
    num := Numerator(omega/BY[#BY]);
    Rnum := Parent(num);
    row := [ MonomialCoefficient(num, Rnum.2^i) : i in [0..(genus - 1)] ];
    Append(~rows, row);
end for;
T := Matrix(rows); TCC := EmbedMatrixExtra(T);

/*
At this point, the period matrix for X is that of Y (taken wrt x^i dx / y), multiplied with T. So PX = T PY. Therefore if A PY = PY R, then
    A T^(-1) T PY = PY R
    T A T^(-1) T PY = T PY R
    T A T^(-1) PX = PX R
and T A T^(-1) is the new tangent representation.
*/

/* Hyperelliptic curve and its period matrix */
Y := HyperellipticCurve(h21(g)); PY := PeriodMatrix(Y);
PX := ChangeRing(TCC, BaseRing(Parent(PY))) * PY;
return PX;

end function;
