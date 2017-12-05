AttachSpec("../../spec");
SetVerbose("EndoCheck", 3);

R<t> := PolynomialRing(Rationals());
F<r> := NumberField(t^2 - t - 1);
R<x> := PolynomialRing(F);
p := 5*x^6 + 10*x^3 - 4*x + 1;
X := HyperellipticCurve(p);
P0 := X ! [0, 1];
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

M := Matrix(F, [
[ r, 0     ],
[ 0, 1 - r ]
]);
M := -M;

print "Field:";
print F;
print "Curve:";
print X;
print "Tangent representation:";
print M;

print "Calculating divisor:";
time test, D := DivisorFromMatrixSplit(X, P0, X, P0, M : LowerBound := 6);
eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
print "Divisor:";
print D;

exit;
