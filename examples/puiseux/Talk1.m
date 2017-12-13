AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 4);

F := Rationals();
R<x> := PolynomialRing(F);
f := -x^5;
h := x^3 + x + 1;
X := HyperellipticCurve(f, h);
P0 := X ! [0, 0, 1];

T := Matrix(F, [
[ -1, -1],
[ -1,  0]
]);

print "Curve:";
print X;
print "Calculating divisor:";
time test, D := DivisorFromMatrixSplit(X, P0, X, P0, T : LowerBound := 6);
eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
print "Divisor:";
print D;

SetVerbose("EndoCheck", 3);
print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixSplit(X, P0, X, P0, T : LowerBound := 11);
R<x,y> := Parent(fs[1]);
print "Cantor representation:";
print fs;

exit;
