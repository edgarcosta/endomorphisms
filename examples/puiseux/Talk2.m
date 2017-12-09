SetVerbose("EndoCheck", 1);

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
time test, D := DivisorFromMatrixSplit(X, P0, X, P0, T : LowerBound := 10);
eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
print "Divisor:";
print D;

print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixSplit(X, P0, X, P0, T : LowerBound := 20);
R<x,y> := Parent(fs[1]);
print "Cantor representation:";
print fs;

exit;
