AttachSpec("../../spec");
SetVerbose("EndoCheck", 3);

R<t> := PolynomialRing(Rationals());
F<r> := NumberField(t^2 + 3);
R<x> := PolynomialRing(F);
p := 24*x^5 + 36*x^4 - 4*x^3 - 12*x^2 + 1;
X := HyperellipticCurve(p);
P0 := X ! [0, 1];
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

M := Matrix(F, [
[ -r, 2*r ],
[  r,   r ]
]);

print "Field:";
print F;
print "Curve:";
print X;
print "Tangent representation:";
print M;

print "Calculating divisor:";
time test, D := DivisorFromMatrixSplit(X, P0, X, P0, M : LowerBound := 28);
eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
print "Divisor:";
print D;

print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixSplit(X, P0, X, P0, M : LowerBound := 104);
R<x,y> := Parent(fs[1]);
print fs;

exit;
