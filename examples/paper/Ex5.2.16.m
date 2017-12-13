/***
 *  An example from the paper
 */

AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 4);

F := FiniteField(4001);
R<x> := PolynomialRing(F);
r := F ! 2845;
i := F ! 3102;

p := x^5 - x^4 + 4*x^3 - 8*x^2 + 5*x - 1;
X := HyperellipticCurve(p);
P0 := X ! [0, i];
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

M := Matrix(F, [
[ 0, r ],
[ r, 0 ]
]);

print "Field:";
print F;
print "Curve:";
print X;
print "Tangent representation:";
print M;

print "Calculating divisor:";
time test, D := DivisorFromMatrix(X, P0, X, P0, M : LowerBound := 8, DivPP1 := true);
eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
print "Divisor:";
print D;
print "Bidegree:";
print Bidegree(X, X, D);

eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
S<x2,x1> := PolynomialRing(F, 2);
res := hom<R -> S | [0, 0, x2, x1]>;
I := ideal<S | res(eqs[#eqs])>;
G := GroebnerBasis(I);

A := AffineSpace(S);
D := Scheme(A, G);
Is := IrreducibleComponents(D);
print "Divisor on PP^1 x PP^1:";
print Is;

print "Calculating Cantor representation...";
time test, fs := CantorFromMatrix(X, P0, X, P0, M : LowerBound := 24);
R<x,y> := Parent(fs[1]);
print fs;

exit;
