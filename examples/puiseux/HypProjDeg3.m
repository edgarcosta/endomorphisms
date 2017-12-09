SetVerbose("EndoCheck", 3);

R<t> := PolynomialRing(Rationals());
F<r> := NumberField(t^2 + 3);
R<x> := PolynomialRing(F);
fX := 60*x^5 + 200*x^4 + 220*x^3 + 89*x^2 + 12*x;
hX := 0;
X := HyperellipticCurve(fX, hX);
P0 := X ! [-1, r, 1];
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

fY := x^3 - 10800*x - 442800;
//fY := x^3 - 6481/3*x - 1042958/27;
Y := HyperellipticCurve(fY);
Q0 := Y ! [1, 0, 0];

M := Matrix(F, [
[2/3, 1]
]);
//M := Matrix(F, [
//[1, 0]
//]);

print "Field:";
print F;
print "Curve X:";
print X;
print "Curve Y:";
print Y;
print "Tangent representation:";
print M;

print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixSplit(X, P0, Y, Q0, M : LowerBound := 24);
R<x,y> := Parent(fs[1]);
print fs;

print "Verification:", CorrespondenceVerifyG1(X, P0, Y, Q0, M, fs);

exit;
