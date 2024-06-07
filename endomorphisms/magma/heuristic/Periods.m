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


/*
function RiemannSurfaceMN(gCC, e : Precision := 30)
return SE_Curve(gCC, e : Prec := Precision);
end function;
*/


function TransformForm(f, T : co := true, contra := false)
    R := Parent(f);
    vars := Matrix([ [ mon ] : mon in MonomialsOfDegree(R, 1) ]);
    if (not co) or contra then
        T := Transpose(T)^(-1);
    end if;
    return Evaluate(f, Eltseq(ChangeRing(T, R) * vars));
end function;


function RandomInvertibleMatrix(n, B)
  D := [ -B..B ];
  repeat
    T := Matrix(Rationals(), n, n, [ Random(D) : i in [1..n^2] ]);
  until Determinant(T) eq 1;
  return T;
end function;

function PeriodMatrixRetry(X, CC: sigma:=false)
  P2 := AmbientSpace(X);
  K := BaseRing(P2);
  assert Dimension(P2) eq 2 and IsProjective(P2); // we expected P2 for the ambient sapce
  if sigma cmpeq false then
    assert Type(K) eq FldRat;
    K := RationalsAsNumberField();
	  X := ChangeRing(X, K);
    P2 := AmbientSpace(X);
	  sigma := InfinitePlaces(K)[1];
  end if;
  F := DefiningPolynomial(X);
  f := DefiningPolynomial(AffinePatch(X, 1));

  T := IdentityMatrix(Rationals(), 3);
  i := 0;
  while true do
    i +:= 1;
    vprintf EndoFind, 3 : "Trying to compute BigPeriodMatrix, attempt = %o\n", i;
    vtime EndoFind, 3:
    try
      RS := RiemannSurface(f, sigma : Precision := Precision(CC));
      break;
    catch e
      vprintf EndoFind, 3 : "Failed %o\n", e;
    end try;
    T := ChangeRing(RandomInvertibleMatrix(3, 2), K);
    newX := PlaneCurve(TransformForm(F, T));
    f := DefiningPolynomial(AffinePatch(newX, 1));
  end while;
  BPM := BigPeriodMatrix(RS);
  P := ChangeRing(BPM, CC);
  // is there a point to keep track of the change of basis on omega_X
  // for genus > 3?
  // if, so why not just look at RS`DFF ?
  if Nrows(P) eq 3 then
    P := ChangeRing(T, BaseRing(P)) * P;
  end if;
  return P, RS;
end function;

intrinsic PeriodMatrix(X::Crv : prec:=false) -> ModMatFldElt
{Returns the period matrix of X. The optional parameter prec is only used if the curve is not given over RationalsExtra/NumberFieldExtra}


F := BaseRing(X);
if not assigned F`iota then
    // converts a curve to be over a NumberFieldExtra if it is not given in that way
    X := CurveExtra(X : prec:=prec);
else
    require prec cmpeq false: "The optional parameter prec can only used if the curve is not given over RationalsExtra/NumberFieldExtra";
end if;
assert assigned F`iota;
CC := Parent(F`iota);
vprint EndoFind : "";
vprint EndoFind : "Calculating period matrix...";
if assigned X`period_matrix and Precision(BaseRing(X`period_matrix)) ge Precision(CC) then
    vprint EndoFind : "using stored period matrix.";
    return ChangeRing(X`period_matrix, CC);
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
    //RS := RiemannSurfaceMN(gCC, 2 : Precision := Precision(CC));
    RS := RiemannSurface(gCC, 2 : Precision := Precision(CC));
    P := ChangeRing(RS`BigPeriodMatrix, CC) / 2;
end if;

if CurveType(X) in [ "plane", "gen" ] then

    /* Determine correct embedding */
    vprintf EndoFind, 3 : "Determining correct embedding...";
    vtime EndoFind, 3:
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
    else
      sigma := false;
      cc := false;
    end if;

    vprintf EndoFind, 3 : "Computing plane model...";
    vtime EndoFind, 3:
    Y := PlaneModel(X);

    P, RS := PeriodMatrixRetry(Y, CC : sigma:=sigma);
    /* Conjugate if needed */
    if cc then
      P := Matrix(CC, [ [ ComplexConjugate(c) : c in Eltseq(row) ] : row in Rows(P) ]);
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
