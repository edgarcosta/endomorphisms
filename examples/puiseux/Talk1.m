AttachSpec("../../spec");
SetVerbose("EndoCheck", 1);

F := Rationals();
R<x> := PolynomialRing(F);
f := -x^5;
h := x^3 + x + 1;
X := HyperellipticCurve(f, h);
P0 := X ! [0, 0];

M := Matrix(F, [
[ -1, -1],
[ -1,  0]
]);

print "Field:";
print F;
print "Curve:";
print X;
print "Tangent representation:";
print M;

print "Calculating divisor:";
time test, D := DivisorFromMatrixSplit(X, P0, X, P0, M : LowerBound := 1, DivPP1 := false);
eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
print "Divisor:";
print D;

exit;

print "Groebner basis for defining equations:";
print GroebnerBasis(ideal<R | eqs>);

print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixSplit(X, P0, X, P0, M : LowerBound := 1);
R<x,y> := Parent(fs[1]);
print "Cantor representation:";
print fs;

exit;
