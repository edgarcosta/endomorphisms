AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 2);

R<t> := PolynomialRing(Rationals());
F<r> := NumberField(t^2 - 2);
R<x> := PolynomialRing(F);
f := x^5 - x^4 + 4*x^3 - 8*x^2 + 5*x - 1;
f := -f;
X := HyperellipticCurve(f);
P0 := X ! [0, 1, 1];

T := Matrix(F, [
[ 1, -r],
[ -r, 1]
]);
//T := -T;

print "Field:";
print F;
print "Curve:";
print X;
print "Tangent representation:";
print T;

print "Calculating divisor:";
time test, D := DivisorFromMatrixAmbientSplit(X, P0, X, P0, T : LowerBound := 1, DivPP1 := false);
eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
print "Divisor:";
print D;

/*
S<x2,x1> := PolynomialRing(F, 2);
h := hom< R -> S | [ 0, 0, x2, x1 ] >;
eqs := DefiningEquations(D)[5..8];
I := EliminationIdeal(Ideal(D), { 3, 4 });
Y := Scheme(AffineSpace(S), [ h(c) : c in GroebnerBasis(I) ]);
Is := [ I : I in IrreducibleComponents(Y) | Dimension(I) eq 1 ];
print Is[1];
*/

print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixAmbientSplit(X, P0, X, P0, T : LowerBound := 36);
R<x,y> := Parent(fs[1]);
print "Cantor representation:";
print fs;

exit;
