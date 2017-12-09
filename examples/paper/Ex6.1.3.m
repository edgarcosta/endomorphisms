/***
 *  An example from the paper
 */

SetVerbose("EndoCheck", 3);

R<t> := PolynomialRing(Rationals());
F<r> := NumberField(t^2 - 2);
R<x> := PolynomialRing(F);
p := x^5 - x^4 + 4*x^3 - 8*x^2 + 5*x - 1;
p *:= -1;
//p *:= -2;
X := HyperellipticCurve(p);
//P0 := X ! [0, r];
P0 := X ! [0, 1];
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
time test, D := DivisorFromMatrixSplit(X, P0, X, P0, M : LowerBound := 8, DivPP1 := true);
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
time test, fs := CantorFromMatrixSplit(X, P0, X, P0, M : LowerBound := 24);
R<x,y> := Parent(fs[1]);
print fs;

exit;
