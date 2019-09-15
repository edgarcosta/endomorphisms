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
forward IsSuperellipticEquation;
forward SuperellipticCompatibility;
forward PeriodMatrixGH;


intrinsic PeriodMatrix(eqsCC::SeqEnum, eqsK::SeqEnum) -> ModMatFldElt, .
{Returns the period matrix of the curve defined by the complex polynomials
eqsCC.}

/* At this point CC has some extra precision */
RCC := Parent(eqsCC[1]); CC := BaseRing(RCC);
if #GeneratorsSequence(RCC) eq 1 then
    if #eqsCC eq 2 then
        fCC, hCC := Explode(eqsCC);
        gCC := (4*fCC + hCC^2) / 4;
    else
        gCC := Explode(eqsCC);
    end if;
    /* We divide by 2 because we integrate with respect to the canonical
     * differential x^i dx / 2y
     * (MN use x^i dx) */
    X := SE_Curve(gCC, 2 : Prec := Precision(CC));
    return ChangeRing(X`BigPeriodMatrix, CC) / 2, X;

elif #GeneratorsSequence(RCC) eq 3 then
    test, fCC, e := IsSuperellipticEquation(eqsCC);
    if test then
        X := SE_Curve(fCC, e : Prec := Precision(CC));
        P := X`BigPeriodMatrix;
        P := SuperellipticCompatibility(P, e);
        return ChangeRing(P, CC), X;
    else
        /* Note: only polynomials over QQ for now */
        F := Explode(eqsK);
        X := PlaneCurve(F); f := DefiningEquation(AffinePatch(X, 1));

        S := RiemannSurface(f : Precision := Precision(CC));
        P := ChangeRing(BigPeriodMatrix(S), CC);
        return P, S;
    end if;

else
    error "No functionality for general curves available";
end if;
end intrinsic;


intrinsic PeriodMatrix(X::SECurve) -> ModMatFldElt
{Returns the period matrix of X.}
F := BaseRing(Parent(X`DefiningPolynomial)); CC := Parent(F`iota);
P := ChangeRing(X`BigPeriodMatrix, CC);

/* Compatibility with Birkenhake--Lange */
g := #Rows(P);
P1 := Submatrix(P, 1,1,   g,g);
P2 := Submatrix(P, 1,g+1, g,g);
P := HorizontalJoin(P2, P1);

return P, X;
end intrinsic;


/* TODO: Also implement this right for generalized hyperelliptic curves for
 * isogeny purposes (right now this is not needed because of the conjugation
 * involved) */
intrinsic PeriodMatrix(X::Crv) -> ModMatFldElt
{Returns the period matrix of X.}

F := BaseRing(X); CC := Parent(F`iota);
vprint EndoFind : "";
vprint EndoFind : "Calculating period matrix...";
if assigned X`period_matrix then
    vprint EndoFind : "using stored period matrix.";
    return X`period_matrix;
end if;

if assigned X`ghpols then
    q, f := Explode(X`ghpols);
    X`period_matrix := PeriodMatrixGH(q, f);
    return X`period_matrix;
end if;

Y := PlaneModel(X);
eqsCC := EmbedCurveEquations(Y); eqsF := DefiningEquations(Y);
P, RS := PeriodMatrix(eqsCC, eqsF);

/* Compatibility with Birkenhake--Lange */
g := #Rows(P);
P1 := Submatrix(P, 1,1,   g,g);
P2 := Submatrix(P, 1,g+1, g,g);
P := HorizontalJoin(P2, P1);

X`period_matrix := ChangeRing(P, CC);
X`riesrf := RS;
vprint EndoFind : "done calculating period matrix.";
return X`period_matrix;

end intrinsic;


/* TODO: Next functions should go since this realization is up to the user once
 * there is a dedicated SE class */
function IsSuperellipticEquation(eqs)
// Returns whether the plane curve defined by eqs is of the form y^e z^* = f
// (x, z). If so, return the inhomogenous form of f along with e.

R<x,y,z> := Parent(eqs[1]);
if #GeneratorsSequence(R) eq 1 then
    return false, 0, 1;
end if;

F := Explode(eqs);
mons := Monomials(F);
monsy := [ mon : mon in mons | Exponents(mon)[2] ne 0 ];
monsxz := [ mon : mon in mons | Exponents(mon)[2] eq 0 ];
if #monsy ne 1 or not &and[ Exponents(mon)[1] eq 0 : mon in monsy ] then
    return false, 0, 1;
end if;

e := Exponents(monsy[1])[2];
S<t> := PolynomialRing(BaseRing(R));
f := &+[ MonomialCoefficient(F, mon) * t^(Exponents(mon)[1]) : mon in monsxz ];
C := MonomialCoefficient(F, monsy[1]);
f := -f/C;
return true, f, e;

end function;


function SuperellipticCompatibility(P, e)
// Transforms the differentials on a superelliptic curve to compensate for
// conventions.
// TODO: Generalize this to apply beyond genus 3. This is a matter of fixing a
// base of differentials. But actually superelliptic curves should be treated
// as a class of their own. Not now: use base provided by Christian.

rowsP := Rows(P);
if #rowsP eq 3 then
    return Matrix([ rowsP[3], rowsP[1], rowsP[2] ]);
end if;
error "Need g = 3 for now";

end function;


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
    p := Evaluate(q, [ R.1, 0, 1 ]); K := NumberFieldExtra(p);
    S := PolynomialRing(K, 3); h := hom< Parent(q) -> S | [ S.1, S.2, S.3 ] >;
    q := h(q); f := h(f);

    S3 := S; KS3 := FieldOfFractions(S3); PP2 := ProjectiveSpace(S3);
    S2 := PolynomialRing(K, 2); KS2 := FieldOfFractions(S2); PP1 := ProjectiveSpace(S2);
    R := PolynomialRing(K); h21 := hom< S2 -> R | [ R.1, 1 ] >;

    Q := Conic(PP2, q); test, P := HasRationalPoint(Q);
end if;
phi := Parametrization(Q, P);
DE := DefiningEquations(phi);
h := hom< Parent(DE[1]) -> S2 | [ S2.1, S2.2 ] >;
DE := [ h(c) : c in DE ];

/* Create hyperelliptic curve */
// TODO: This could give incompatibilities! Not yet though.
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
